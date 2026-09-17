#define NS_LOG_APPEND_CONTEXT                                                                       \
    do                                                                                              \
    {                                                                                               \
        std::clog << "[ CellId " << GetCellId() << ", bwpId " << GetBwpId() << "] ";                \
    } while (false);                                                    

#include "nr-mac-scheduler-tdma-sdma.h"
#include "nr-mac-scheduler-ue-info-rr.h"

#include "ns3/log.h"

#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <complex>
#include <limits>

#include "ns3/boolean.h"
#include "ns3/enum.h"
#include "ns3/integer.h"
#include "ns3/log.h"
#include "ns3/uinteger.h"
#include "ns3/double.h"

#include <Eigen/Dense>


namespace ns3
{

NS_LOG_COMPONENT_DEFINE("NrMacSchedulerTdmaSdma");
NS_OBJECT_ENSURE_REGISTERED(NrMacSchedulerTdmaSdma);

using UePtr = std::shared_ptr<NrMacSchedulerUeInfo>;

// ---------- helpers ----------------------------------------------------------
namespace
{
using Cd  = std::complex<double>;
using Mat = Eigen::MatrixXcd;
using Vec = Eigen::VectorXcd;

// You already have this; included for completeness. Pulls one subband page
// out of the ComplexMatrixArray into an Nrx x Ntx Eigen matrix.
inline Mat
PageToEigen(const Ptr<const ComplexMatrixArray>& H, size_t page)
{
    const size_t nrx = H->GetNumRows();
    const size_t ntx = H->GetNumCols();
    Mat M(nrx, ntx);
    for (size_t r = 0; r < nrx; ++r)
        for (size_t c = 0; c < ntx; ++c)
            M(r, c) = (*H)(r, c, page);
    return M;
}

// Per-UE wideband TRANSMIT covariance  R = sum_p H_p^H H_p   (Ntx x Ntx).
// Hermitian PSD. Carries Sigma^2 and the tx singular vectors V (R = V Sigma^2 V^H),
// which is exactly what both grouping metrics consume. Summing covariances
// (not channels) avoids per-page phase cancellation.
inline Mat
TxCovariance(const Ptr<const ComplexMatrixArray>& H)
{
    const size_t ntx   = H->GetNumCols();
    const size_t pages = H->GetNumPages();
    Mat R = Mat::Zero(ntx, ntx);
    for (size_t p = 0; p < pages; ++p)
    {
        const Mat Hp = PageToEigen(H, p);   // Nrx x Ntx
        R.noalias() += Hp.adjoint() * Hp;   // Ntx x Ntx
    }
    return R;
}

// Dominant tx direction = principal eigenvector of R (unit norm),
// and its eigenvalue (= dominant energy, used to pick group seeds).
inline void
DominantTxDir(const Mat& R, Vec& dirOut, double& lambdaOut)
{
    Eigen::SelfAdjointEigenSolver<Mat> es(R);          // ascending eigenvalues
    const Eigen::Index last = R.cols() - 1;
    lambdaOut = es.eigenvalues()(last);
    dirOut    = es.eigenvectors().col(last);           // already unit norm
}

// Cached per-UE grouping descriptor.
struct UeGroupInfo
{
    Mat    R;        // tx covariance (Ntx x Ntx)
    Vec    dir;      // dominant tx direction (unit, Ntx x 1)
    double lambda;   // dominant energy (seed ordering)
    double trR2;     // trace(R^2) = ||R||_F^2  (precomputed for collinearity)
};

// ---- (iii) Subspace collinearity in [0,1] -------------------------- //
// c = tr(Ri Rj) / sqrt( tr(Ri^2) tr(Rj^2) ).
// 0 => tx subspaces orthogonal (ideal to co-schedule); 1 => identical.
// trace(Ri Rj) is real for Hermitian Ri,Rj.
inline double
SubspaceCollinearity(const UeGroupInfo& a, const UeGroupInfo& b)
{
    const double num = (a.R * b.R).trace().real();
    const double den = std::sqrt(a.trR2 * b.trR2);
    if (den <= 0.0) return 1.0;
    double c = num / den;
    return std::clamp(c, 0.0, 1.0);
}

// ---- (ii) Matrix volume ratio (exact, Gram-Schmidt) ---------------- //
// rho = || g - Q Q^H g || / ||g||  =  sin(angle between g and span(Q)).
// Q has orthonormal columns spanning the current group's tx directions.
// rho in [0,1]: 1 = orthogonal to the group (adds full volume), 0 = collinear.
inline double
VolumeRatio(const Mat& Q, const Vec& g)
{
    const double gn = g.norm();
    if (gn <= 0.0) return 0.0;
    if (Q.cols() == 0) return 1.0;                 // empty group: full volume
    const Vec resid = g - Q * (Q.adjoint() * g);   // orthogonal component
    return std::clamp(resid.norm() / gn, 0.0, 1.0);
}

// ---- (ii) AM-GM cheap approximation of the volume ratio ------------ //
// Approximates the determinant/volume gain by the mean squared correlation
// to existing members instead of the exact orthogonal residual:
//   rho_approx = sqrt( 1 - mean_m |g^H g_m|^2 )   (g, g_m unit).
// Looser but O(K) and no orthogonalization (the [144]/[262]-style proxy).
inline double
VolumeRatioAmGm(const std::vector<Vec>& groupDirs, const Vec& g)
{
    if (groupDirs.empty()) return 1.0;
    double acc = 0.0;
    for (const Vec& gm : groupDirs)
    {
        const Cd ip = g.dot(gm);          // g^H g_m (Eigen dot conjugates lhs)
        acc += std::norm(ip);             // |.|^2
    }
    acc /= static_cast<double>(groupDirs.size());
    return std::sqrt(std::max(0.0, 1.0 - acc));
}

// Append a (re-orthonormalized) direction as a new column of Q.
inline void
AppendOrthoColumn(Mat& Q, const Vec& g)
{
    Vec v = g;
    if (Q.cols() > 0)
        v -= Q * (Q.adjoint() * g);       // remove component already in span
    const double n = v.norm();
    if (n <= 1e-12) return;               // numerically dependent: skip
    v /= n;
    Mat Qn(Q.rows(), Q.cols() + 1);
    if (Q.cols() > 0) Qn.leftCols(Q.cols()) = Q;
    Qn.col(Q.cols()) = v;
    Q = std::move(Qn);
}

// Per-UE result shared by both precoders, consumed by FinalizeAllocation().
struct UeAlloc
{
    std::shared_ptr<NrMacSchedulerUeInfo> ue;
    Ptr<const ComplexMatrixArray> H;          // [Nrx x Ntx x Pages]
    uint32_t rank = 0;
    std::vector<Eigen::MatrixXcd> Wp;          // per page: Ntx  x rank (precoder)
    std::vector<Eigen::MatrixXcd> Up;          // per page: Nrx  x rank (rx combiner)
};

// Sum_p W_p W_p^H over RBs. For a rank-1 precoder (W_p = w_p, Ntx x 1) this is
// Sum_f w_f w_f^H — exactly the phase-invariant direction covariance.
inline Mat PrecoderCov(const Ptr<const ComplexMatrixArray>& W)
{
    const size_t nTx = W->GetNumRows();
    const size_t nL  = W->GetNumCols();   // == 1 for rank-1
    const size_t nRb = W->GetNumPages();
    Mat R = Mat::Zero(nTx, nTx);
    for (size_t p = 0; p < nRb; ++p)
        for (size_t l = 0; l < nL; ++l)
        {
            Vec w(nTx);
            for (size_t r = 0; r < nTx; ++r)
                w(r) = (*W)(r, l, p);     // adapt ()-indexing to your ComplexMatrixArray
            R += w * w.adjoint();
        }
    return R;
}

// Build the per-UE cache. Free template so a const member can supply channel
// access via a lambda; UeGroupInfo never escapes this .cc.
template <class GetChan>
inline std::map<uint16_t, UeGroupInfo>
BuildGroupingCache(const std::vector<UePtr>& ues, GetChan getChan)
{
    std::map<uint16_t, UeGroupInfo> cache;
    for (const auto& ue : ues)
    {
        if (!ue) continue;
        Ptr<const ComplexMatrixArray> H = getChan(ue);
        if (!H || H->GetNumPages() == 0) continue;
        UeGroupInfo gi;
        gi.R    = TxCovariance(H);
        // gi.R    = PrecoderCov(ue->m_dlPrecMats);
        gi.trR2 = gi.R.squaredNorm();
        DominantTxDir(gi.R, gi.dir, gi.lambda);
        cache.emplace(ue->m_rnti, std::move(gi));
    }
    return cache;
}
 
// Order UEs by dominant energy (strongest user seeds each group).
inline std::vector<UePtr>
SortByEnergy(const std::vector<UePtr>& ues,
             const std::map<uint16_t, UeGroupInfo>& cache)
{
    std::vector<UePtr> out = ues;
    std::stable_sort(out.begin(), out.end(),
        [&](const UePtr& a, const UePtr& b) {
            auto ia = cache.find(a->m_rnti);
            auto ib = cache.find(b->m_rnti);
            double la = (ia != cache.end()) ? ia->second.lambda : 0.0;
            double lb = (ib != cache.end()) ? ib->second.lambda : 0.0;
            return la > lb;
        });
    return out;
}

// ============================================================================
//  K-MEANS user clustering on dominant tx directions
//  ---------------------------------------------------------------------------
//  Operates on the complex projective space CP^{Ntx-1} (each UE = a 1-D beam).
//  Distance : squared chordal d^2(g,c) = 1 - |g^H c|^2  (g,c unit norm) in [0,1]
//             0 = collinear beams, 1 = orthogonal. Phase invariant.
//  Centroid : principal eigenvector of S = sum_k g_k g_k^H  (the Frechet mean
//             of a set of 1-D subspaces under the chordal metric).
//  Init     : deterministic farthest-first (k-means++ without randomness) so
//             simulations stay reproducible run-to-run.
// ============================================================================

// Squared chordal distance between two unit vectors (phase invariant).
inline double
ChordalDistSq(const Vec& g, const Vec& c)
{
    const Cd ip = c.dot(g);                       // c^H g (Eigen dot conjugates lhs)
    return std::max(0.0, 1.0 - std::norm(ip));    // 1 - |g^H c|^2
}

// Cluster centroid = principal eigenvector of sum_{k in members} g_k g_k^H.
inline Vec
ClusterCentroid(const std::vector<Vec>& dirs, const std::vector<int>& members)
{
    const Eigen::Index n = dirs.front().rows();
    Mat S = Mat::Zero(n, n);
    for (int idx : members)
        S.noalias() += dirs[idx] * dirs[idx].adjoint();
    Eigen::SelfAdjointEigenSolver<Mat> es(S);     // ascending eigenvalues
    return es.eigenvectors().col(n - 1);          // principal eigenvector (unit norm)
}

// Deterministic farthest-first seeding (k-means++ without an RNG).
// First seed = dirs[0] (caller passes them in dominant-energy order),
// each subsequent seed = the UE maximizing its min chordal distance to the
// already chosen seeds.
inline std::vector<Vec>
KMeansFarthestInit(const std::vector<Vec>& dirs, uint32_t K)
{
    std::vector<Vec> centroids;
    centroids.push_back(dirs.front());                       // strongest UE first
    std::vector<double> dmin(dirs.size(), std::numeric_limits<double>::max());
    while (centroids.size() < K)
    {
        const Vec& last = centroids.back();
        double bestVal = -1.0;
        size_t bestIdx = 0;
        for (size_t i = 0; i < dirs.size(); ++i)
        {
            dmin[i] = std::min(dmin[i], ChordalDistSq(dirs[i], last));
            if (dmin[i] > bestVal) { bestVal = dmin[i]; bestIdx = i; }
        }
        if (bestVal <= 1e-12) break;                         // remaining points coincide
        centroids.push_back(dirs[bestIdx]);
    }
    return centroids;
}

// Lloyd's algorithm on CP^{Ntx-1}. Returns assignment[i] in [0,K).
inline std::vector<int>
KMeansClusterDirections(const std::vector<Vec>& dirs, uint32_t K, uint32_t maxIter)
{
    const size_t N = dirs.size();
    std::vector<int> assign(N, 0);
    if (N == 0) return assign;
    K = std::min<uint32_t>(K, static_cast<uint32_t>(N));
    if (K <= 1) return assign;                               // single cluster

    std::vector<Vec> centroids = KMeansFarthestInit(dirs, K);
    K = static_cast<uint32_t>(centroids.size());

    for (uint32_t it = 0; it < maxIter; ++it)
    {
        // ---- assignment step ----
        bool changed = false;
        for (size_t i = 0; i < N; ++i)
        {
            double best = std::numeric_limits<double>::max();
            int    bestK = 0;
            for (uint32_t k = 0; k < K; ++k)
            {
                const double d = ChordalDistSq(dirs[i], centroids[k]);
                if (d < best) { best = d; bestK = static_cast<int>(k); }
            }
            if (assign[i] != bestK) { assign[i] = bestK; changed = true; }
        }
        if (!changed) break;                                 // converged

        // ---- update step ----
        std::vector<std::vector<int>> members(K);
        for (size_t i = 0; i < N; ++i)
            members[assign[i]].push_back(static_cast<int>(i));
        for (uint32_t k = 0; k < K; ++k)
            if (!members[k].empty())
                centroids[k] = ClusterCentroid(dirs, members[k]);
        // empty clusters keep their previous centroid
    }
    return assign;
}

} // namespace

TypeId
NrMacSchedulerTdmaSdma::GetTypeId()
{
    static TypeId tid = TypeId("ns3::NrMacSchedulerTdmaSdma")
                        .SetParent<NrMacSchedulerNs3>()
                        .AddAttribute("maxSdmaUsers",
                                        "Maximum number of allowed users to be grouped in SDMA, default = 4 (5G-LENA has the limit of max 4 streams per transmission)",
                                        UintegerValue(4),
                                        MakeUintegerAccessor(&NrMacSchedulerTdmaSdma::m_maxSdmaUsers),
                                        MakeUintegerChecker<uint32_t>(1, 4))
                        .AddAttribute("compatThres",
                                        "Threshold for compatibility of two users based on their precoding matrix. Default = 0.001",
                                        DoubleValue(3e-1),
                                        MakeDoubleAccessor(&NrMacSchedulerTdmaSdma::m_compatThresh),
                                        MakeDoubleChecker<double>(0.0, 1.0))
                        .AddAttribute("FullCsiBaseline", "Enable full-CSI baseline path.",
                                        BooleanValue(false), MakeBooleanAccessor(&NrMacSchedulerTdmaSdma::m_fullCsiBaseline), MakeBooleanChecker())
                        .AddAttribute("AggregateSvd", "false=ZF, true=aggregate-SVD joint precoder.",
                                        BooleanValue(false), MakeBooleanAccessor(&NrMacSchedulerTdmaSdma::m_aggregateSvd), MakeBooleanChecker())
                        .AddAttribute("MaxStreams", "Max streams/symbol; 0 = unconstrained (Part 1), 4 = Part 2.",
                                        UintegerValue(4), MakeUintegerAccessor(&NrMacSchedulerTdmaSdma::m_maxStreams), MakeUintegerChecker<uint32_t>(0,64))
                        .AddAttribute("MaxRankPerUe", "Max ZF streams per UE in the full-CSI baseline.",
                                        UintegerValue(1), MakeUintegerAccessor(&NrMacSchedulerTdmaSdma::m_maxRankPerUe),
                                        MakeUintegerChecker<uint32_t>(1, 4))
                        .AddAttribute("SinrBackoffDb", "SINR backoff (dB) applied before MCS selection (BLER margin).",
                                        DoubleValue(3.0), MakeDoubleAccessor(&NrMacSchedulerTdmaSdma::m_sinrBackoffDb),
                                        MakeDoubleChecker<double>(0.0, 30.0))
                        .AddAttribute("DlNoisePowerLin", "Linear noise power reference for MCS recomputation.",
                                        DoubleValue(1.0), MakeDoubleAccessor(&NrMacSchedulerTdmaSdma::m_dlNoisePowerLin), MakeDoubleChecker<double>(0.0))
                        .AddAttribute("ForceDlMcs", "Debug: force this DL MCS in the baseline (-1 = off).",
                                        IntegerValue(-1), MakeIntegerAccessor(&NrMacSchedulerTdmaSdma::m_forceDlMcs),
                                        MakeIntegerChecker<int>(-1, 28))
                        .AddAttribute("GroupingMetric",
                                    "SDMA user-grouping criterion: 0=ZF-dir, 1=collinearity, 2=volume, 3=K-means",
                                    UintegerValue(0),
                                    MakeUintegerAccessor(&NrMacSchedulerTdmaSdma::m_groupingMetric),
                                    MakeUintegerChecker<uint8_t>(0, 3))
                        .AddAttribute("CollinearityThreshold",
                                    "Max allowed pairwise subspace collinearity to co-group [0,1]",
                                    DoubleValue(0.10),
                                    MakeDoubleAccessor(&NrMacSchedulerTdmaSdma::m_collinearityThr),
                                    MakeDoubleChecker<double>(0.0, 1.0))
                        .AddAttribute("VolumeRatioThreshold",
                                    "Min volume ratio (sin angle) to add a user to a group [0,1]",
                                    DoubleValue(0.50),
                                    MakeDoubleAccessor(&NrMacSchedulerTdmaSdma::m_volumeThr),
                                    MakeDoubleChecker<double>(0.0, 1.0))
                        .AddAttribute("VolumeUseAmGm",
                                    "Use the cheap AM-GM volume approximation instead of exact",
                                    BooleanValue(false),
                                    MakeBooleanAccessor(&NrMacSchedulerTdmaSdma::m_volumeUseAmGm),
                                    MakeBooleanChecker())
                        .AddAttribute("KMeansClusters",
                                    "Number of K-means clusters (0 = use maxSdmaUsers / co-scheduling order)",
                                    UintegerValue(0),
                                    MakeUintegerAccessor(&NrMacSchedulerTdmaSdma::m_kmeansClusters),
                                    MakeUintegerChecker<uint32_t>(0, 64))
                        .AddAttribute("KMeansMaxIter",
                                    "Maximum Lloyd iterations for K-means user clustering",
                                    UintegerValue(20),
                                    MakeUintegerAccessor(&NrMacSchedulerTdmaSdma::m_kmeansMaxIter),
                                    MakeUintegerChecker<uint32_t>(1, 1000));
    return tid;
    return tid;
}

NrMacSchedulerTdmaSdma::NrMacSchedulerTdmaSdma()
{
}

NrMacSchedulerTdmaSdma::~NrMacSchedulerTdmaSdma()
{

}

std::vector<NrMacSchedulerNs3::UePtrAndBufferReq>
NrMacSchedulerTdmaSdma::GetUeVectorFromActiveUeMap(const NrMacSchedulerNs3::ActiveUeMap& activeUes)
{
    std::vector<UePtrAndBufferReq> ueVector;
    for (const auto& el : activeUes)
    {
        uint64_t size = ueVector.size();
        GetSecond GetUeVector;
        for (const auto& ue : GetUeVector(el))
        {
            ueVector.emplace_back(ue);
        }
        NS_ASSERT(size + GetUeVector(el).size() == ueVector.size());
    }
    return ueVector;
}

bool
NrMacSchedulerTdmaSdma::IsCompatibleDl(const UePtrAndBufferReq& ueA,
                                        const UePtrAndBufferReq& ueB) const
{
    // Get the UE Ptr of the two UEs
    const auto& a = ueA.first;
    const auto& b = ueB.first;

    if (!a || !b)
    {
        NS_LOG_DEBUG("One of the UE pointer's is null!");
        return false;
    }

    // Simple assumption: Only rank 1 grouping
    if (a->m_dlRank != 1 || b->m_dlRank != 1)
    {
        NS_LOG_DEBUG("Only rank 1 grouping is allowed");
        return false;
    }

    if (a->m_dlPrecMats == nullptr || b->m_dlPrecMats == nullptr)
    {
        NS_LOG_DEBUG("One of the precoding matrices is null");
        return false;
    }

    double res = 0.0;
    auto Pages = a->m_dlPrecMats->GetNumPages();

    for (size_t page = 0; page < Pages; ++page)
    {
        auto w_a = a->m_dlPrecMats->ExtractPage(page);
        auto w_b = b->m_dlPrecMats->ExtractPage(page);

        std::complex<double> sa = 0.0, sb = 0.0;
        for (uint32_t r = 0; r < w_a.GetNumRows(); ++r)
        {
            sa += std::conj(w_a(r,0)) * w_a(r,0);
            sb += std::conj(w_b(r,0)) * w_b(r,0);
        }
        double na = std::sqrt(sa.real());
        double nb = std::sqrt(sb.real());

        auto prod= w_a.HermitianTranspose() * w_b;
        std::complex<double> s = prod(0, 0) / (na * nb);

        res += std::norm(s);
    }

    res /= Pages;

    NS_LOG_UNCOND("The result between " << a->m_rnti << " and " << b->m_rnti << " is " << res << " and the threshold is " << m_compatThresh << ", Compatible = " << (res <= m_compatThresh));

    return res <= m_compatThresh;
}

void
NrMacSchedulerTdmaSdma::ComputeDlGroupPrecoding(const std::vector<UePtr>& group) const
{
    NS_LOG_UNCOND("THE VALV IS = " << m_fullCsiBaseline);
    if (m_fullCsiBaseline)
    {
        if (m_aggregateSvd) ComputeAggregateSvdPrecoding(group, m_maxStreams);
        else                ComputeFullCsiSvdPrecoding(group, m_maxStreams);
    }
    else
    {
        ComputeGroupZfPrecoding(group);
    }
}

// ---- (iii) Greedy grouping by subspace collinearity ---------------- //
// Add a user to the current group iff its MAX pairwise collinearity with the
// members stays <= m_collinearityThr. Cap group size at maxStreams.
std::vector<std::vector<UePtr>>
NrMacSchedulerTdmaSdma::GroupByCollinearity(const std::vector<UePtr>& ues,
                                            uint32_t maxStreams) const
{
    auto cache = BuildGroupingCache(ues,
        [this](const UePtr& u) { return GetDlChannelMatrix(u); });
    auto pool  = SortByEnergy(ues, cache);
    std::vector<std::vector<UePtr>> groups;
 
    std::vector<bool> used(pool.size(), false);
    for (size_t s = 0; s < pool.size(); ++s)
    {
        if (used[s] || !cache.count(pool[s]->m_rnti)) continue;
        std::vector<UePtr> group{pool[s]};
        used[s] = true;
 
        for (size_t c = 0; c < pool.size(); ++c)
        {
            if (used[c] || group.size() >= maxStreams) continue;
            auto itc = cache.find(pool[c]->m_rnti);
            if (itc == cache.end()) continue;
 
            double worst = 0.0;
            for (const auto& m : group)
                worst = std::max(worst,
                                 SubspaceCollinearity(cache.at(m->m_rnti), itc->second));
 
            if (worst <= m_collinearityThr)
            {
                group.push_back(pool[c]);
                used[c] = true;
            }
        }
        groups.push_back(std::move(group));
    }
    return groups;
}


std::vector<std::vector<UePtr>>
NrMacSchedulerTdmaSdma::GroupByVolume(const std::vector<UePtr>& ues,
                                      uint32_t maxStreams) const
{
    auto cache = BuildGroupingCache(ues,
        [this](const UePtr& u) { return GetDlChannelMatrix(u); });
    auto pool  = SortByEnergy(ues, cache);
    std::vector<std::vector<UePtr>> groups;
 
    std::vector<bool> used(pool.size(), false);
    for (size_t s = 0; s < pool.size(); ++s)
    {
        if (used[s] || !cache.count(pool[s]->m_rnti)) continue;
 
        std::vector<UePtr> group{pool[s]};
        used[s] = true;
 
        const Vec& seedDir = cache.at(pool[s]->m_rnti).dir;
        Mat Q(seedDir.rows(), 0);
        AppendOrthoColumn(Q, seedDir);
        std::vector<Vec> groupDirs{seedDir};
 
        while (group.size() < maxStreams)
        {
            int    best    = -1;
            double bestRho = m_volumeThr;
            for (size_t c = 0; c < pool.size(); ++c)
            {
                if (used[c]) continue;
                auto itc = cache.find(pool[c]->m_rnti);
                if (itc == cache.end()) continue;
 
                const double rho = m_volumeUseAmGm
                    ? VolumeRatioAmGm(groupDirs, itc->second.dir)
                    : VolumeRatio(Q, itc->second.dir);
 
                if (rho > bestRho) { bestRho = rho; best = static_cast<int>(c); }
            }
            if (best < 0) break;
 
            const Vec& bd = cache.at(pool[best]->m_rnti).dir;
            AppendOrthoColumn(Q, bd);
            groupDirs.push_back(bd);
            group.push_back(pool[best]);
            used[best] = true;
        }
        groups.push_back(std::move(group));
    }
    return groups;
}
 
// ---- (iv) K-means user clustering (cluster-then-pair) -------------- //
// 1. Cluster UEs by spatial similarity (chordal K-means on dominant tx dirs).
//    Users in the same cluster share a beam => would interfere if co-scheduled.
// 2. Build each SDMA group by taking one UE from each cluster, so group members
//    come from DIFFERENT clusters => spatially separated => good to co-schedule.
// K = co-scheduling order (m_maxSdmaUsers by default), capped by maxStreams.
std::vector<std::vector<UePtr>>
NrMacSchedulerTdmaSdma::GroupByKMeans(const std::vector<UePtr>& ues,
                                      uint32_t maxStreams) const
{
    auto cache = BuildGroupingCache(ues,
        [this](const UePtr& u) { return GetDlChannelMatrix(u); });

    // Keep only UEs with a valid descriptor, in dominant-energy order
    // (strong users seed clusters and are paired first).
    auto pool = SortByEnergy(ues, cache);
    std::vector<UePtr> valid;
    std::vector<Vec>   dirs;
    valid.reserve(pool.size());
    dirs.reserve(pool.size());
    for (const auto& ue : pool)
    {
        auto it = cache.find(ue->m_rnti);
        if (it == cache.end()) continue;
        valid.push_back(ue);
        dirs.push_back(it->second.dir);
    }
    if (valid.empty()) return {};

    // Number of clusters = desired co-scheduling order, capped by stream budget.
    const uint32_t desiredK = (m_kmeansClusters == 0) ? m_maxSdmaUsers : m_kmeansClusters;
    uint32_t K = std::min<uint32_t>({desiredK,
                                     std::max<uint32_t>(maxStreams, 1),
                                     static_cast<uint32_t>(valid.size())});
    K = std::max<uint32_t>(K, 1);

    const std::vector<int> assign =
        KMeansClusterDirections(dirs, K, m_kmeansMaxIter);

    // Bucket UEs by cluster, preserving energy order within each cluster.
    std::vector<std::vector<UePtr>> clusters(K);
    for (size_t i = 0; i < valid.size(); ++i)
        clusters[assign[i]].push_back(valid[i]);

    // ---- form SDMA groups: row t = one UE from each cluster ----
    // group_t = { clusters[0][t], clusters[1][t], ... } capped at maxStreams.
    size_t maxClusterSize = 0;
    for (const auto& c : clusters)
        maxClusterSize = std::max(maxClusterSize, c.size());

    std::vector<std::vector<UePtr>> groups;
    std::vector<UePtr> leftovers;          // UEs rejected by the collinearity guard

    for (size_t t = 0; t < maxClusterSize; ++t)
    {
        std::vector<UePtr> group;
        for (uint32_t k = 0; k < K && group.size() < maxStreams; ++k)
            if (t < clusters[k].size())
                group.push_back(clusters[k][t]);

        // Collinearity guard: always keep the seed (group[0]); test each later
        // member against everyone already kept. Reject -> re-home as a singleton,
        // never drop, so the partition stays complete.
        for (size_t i = 1; i < group.size(); )
        {
            bool ok = true;
            for (size_t j = 0; j < i; ++j)
                if (SubspaceCollinearity(cache.at(group[j]->m_rnti),
                                         cache.at(group[i]->m_rnti)) > m_collinearityThr)
                { ok = false; break; }

            if (ok) { ++i; }
            else    { leftovers.push_back(group[i]); group.erase(group.begin() + i); }
        }

        if (!group.empty()) groups.push_back(std::move(group));
    }

    // Re-home every rejected UE as its own TDMA group -> partition is complete.
    for (const auto& ue : leftovers)
        groups.push_back({ue});

    return groups;
}
 
std::vector<std::vector<UePtr>>
NrMacSchedulerTdmaSdma::GroupUesForSdma(const std::vector<UePtr>& ues) const
{
    const uint32_t cap = (m_maxStreams == 0)
                             ? static_cast<uint32_t>(ues.size())
                             : std::min<uint32_t>(m_maxStreams,
                                                  static_cast<uint32_t>(ues.size()));
    switch (m_groupingMetric)
    {
        case GROUP_COLINEAR: return GroupByCollinearity(ues, cap);
        case GROUP_VOLUME:       return GroupByVolume(ues, cap);
        case GROUP_KMEANS:       return GroupByKMeans(ues, cap);
        default:                 return {};   // GROUP_ZF_DIR handled by caller
    }
}


NrMacSchedulerTdmaSdma::BeamSymbolMap
NrMacSchedulerTdmaSdma::AssignRBGTDMA(uint32_t symAvail,
                                  const ActiveUeMap& activeUe,
                                  const std::string& type,
                                  const BeforeSchedFn& BeforeSchedFn,
                                  const GetCompareUeFn& GetCompareFn,
                                  const GetTBSFn& GetTBSFn,
                                  const GetRBGFn& GetRBGFn,
                                  const GetSymFn& GetSymFn,
                                  const AfterSuccessfulAssignmentFn& SuccessfulAssignmentFn,
                                  const AfterUnsuccessfulAssignmentFn& UnSuccessfulAssignmentFn,
                                  const CallNotifyFn& callNotifyFn) const
{
    NS_LOG_FUNCTION(this);

    std::vector<UePtrAndBufferReq> ueVector = GetUeVectorFromActiveUeMap(activeUe);

    uint32_t resources = symAvail;
    FTResources assigned(0, 0);

    const std::vector<bool> notchedRBGsMask = type == "DL" ? GetDlBitmask() : GetUlBitmask();
    uint32_t numOfAssignableRbgs = std::count(notchedRBGsMask.begin(), notchedRBGsMask.end(), true);
    NS_ASSERT(numOfAssignableRbgs > 0);

    std::set<uint32_t> remainingRbgSet;
    for (size_t i = 0; i < notchedRBGsMask.size(); i++)
    {
        if (notchedRBGsMask.at(i))
        {
            remainingRbgSet.emplace(i);
        }
    }

    for (auto& ue : ueVector)
    {
        BeforeSchedFn(ue, FTResources(numOfAssignableRbgs, 1));
    }

    GetFirst GetUe;

    auto NeedsMoreResources = [&](const UePtrAndBufferReq& u) -> bool {
        uint32_t bufQueueSize = u.second;
        return (GetTBSFn(GetUe(u)) >= std::max(bufQueueSize, 10U));
    };

    // >>> -------- baseline mode + ZF orthogonality grouping setup --------------
    const bool zfBaseline = (type == "DL" && m_fullCsiBaseline);           // >>>
    const bool noLimit    = zfBaseline && (m_maxStreams == 0);             // >>>

    // Cache each UE's per-subband dominant channel direction (top right singular
    // vector of H_k). Only needed for the orthogonality test (limited ZF).
    std::unordered_map<uint16_t, Eigen::MatrixXcd> dir;                    // >>>
    uint32_t Ntx = 0;                                                      // >>>
    if (zfBaseline && !noLimit)                                            // >>>
    {                                                                      // >>>
        for (const auto& ue : ueVector)                                    // >>>
        {                                                                  // >>>
            auto H = GetDlChannelMatrix(GetUe(ue));                        // >>>
            if (!H || H->GetNumPages() == 0 || H->GetNumCols() == 0)       // >>>
                continue;                                                  // >>>
            const size_t pg = H->GetNumPages();                            // >>>
            Ntx = H->GetNumCols();                                         // >>>
            Eigen::MatrixXcd G(Ntx, pg);                                   // >>>
            for (size_t p = 0; p < pg; ++p)                                // >>>
            {                                                              // >>>
                Eigen::JacobiSVD<Eigen::MatrixXcd> svd(PageToEigen(H, p), // >>>
                                                       Eigen::ComputeThinV);// >>>
                G.col(p) = svd.matrixV().col(0);                           // >>>
            }                                                              // >>>
            dir[GetUe(ue)->m_rnti] = std::move(G);                         // >>>
        }                                                                  // >>>
    }                                                                     // >>>

    auto IsCompatibleZf = [&](uint16_t a, uint16_t b) -> bool {            // >>>
        auto ia = dir.find(a), ib = dir.find(b);                          // >>>
        if (ia == dir.end() || ib == dir.end()) return false;             // >>>
        const auto& Ga = ia->second; const auto& Gb = ib->second;         // >>>
        const Eigen::Index P = std::min(Ga.cols(), Gb.cols());            // >>>
        double rho = 0.0;                                                 // >>>
        for (Eigen::Index p = 0; p < P; ++p)                              // >>>
            rho += std::norm((Ga.col(p).adjoint() * Gb.col(p))(0, 0));    // >>>
        rho /= std::max<Eigen::Index>(P, 1);                              // >>>
        return rho <= m_compatThresh;   // small = orthogonal = good ZF   // >>>
    };                                                                     // >>>

    // group-size cap by mode
    const uint32_t groupCap =                                             // >>>
          noLimit    ? std::numeric_limits<uint32_t>::max()               // >>>
        : zfBaseline ? (Ntx ? std::min<uint32_t>(m_maxStreams, Ntx)       // >>>
                            : m_maxStreams)                               // >>>
                     : m_maxSdmaUsers;                                    // >>>
    // >>> ----------------------------------------------------------------------

    // ============================================================
    // PHASE 1: build fixed SDMA groups (each UE belongs to exactly one group)
    // ============================================================
    SortUeVector(&ueVector, GetCompareFn);

    std::vector<std::vector<UePtrAndBufferReq>> groups;
    std::unordered_map<uint16_t, size_t> rntiToGroup;
    if (m_groupingMetric == GROUP_ZF_DIR)
    {
        for (auto seedIt = ueVector.begin(); seedIt != ueVector.end(); ++seedIt)
        {
            const uint16_t seedRnti = GetUe(*seedIt)->m_rnti;
            if (rntiToGroup.count(seedRnti) != 0)
            {
                continue;
            }

            std::vector<UePtrAndBufferReq> group;
            group.reserve(ueVector.size());                                   // >>> (was m_maxSdmaUsers)
            group.push_back(*seedIt);
            rntiToGroup[seedRnti] = groups.size();

            for (auto candIt = ueVector.begin();
                candIt != ueVector.end() && group.size() < groupCap;         // >>> (was m_maxSdmaUsers)
                ++candIt)
            {
                const uint16_t candRnti = GetUe(*candIt)->m_rnti;
                if (candRnti == seedRnti || rntiToGroup.count(candRnti) != 0)
                {
                    continue;
                }

                // >>> compatibility test depends on mode
                bool compatibleWithWholeGroup = true;                         // >>>
                if (noLimit)                                                  // >>>
                {                                                             // >>>
                    compatibleWithWholeGroup = true;       // all together    // >>>
                }                                                             // >>>
                else if (zfBaseline)                                          // >>>
                {                                                             // >>>
                    for (const auto& member : group)                         // >>>
                        if (!IsCompatibleZf(GetUe(member)->m_rnti, candRnti)) // >>>
                        { compatibleWithWholeGroup = false; break; }          // >>>
                }                                                             // >>>
                else                                                          // >>>
                {                                                             // >>>
                    for (const auto& member : group)                         // >>>
                        if (!IsCompatibleDl(member, *candIt))                 // >>>
                        { compatibleWithWholeGroup = false; break; }          // >>>
                }                                                             // >>>

                if (compatibleWithWholeGroup)
                {
                    rntiToGroup[candRnti] = groups.size();
                    group.push_back(*candIt);
                }
            }

            groups.push_back(std::move(group));
        }
    }
    else
    {
        std::vector<UePtr> bare;
        bare.reserve(ueVector.size());
        for (const auto& u : ueVector) bare.push_back(GetUe(u));

        auto bareGroups = GroupUesForSdma(bare);

        std::unordered_map<uint16_t, UePtrAndBufferReq> byRnti;
        for (const auto& u : ueVector) byRnti.emplace(GetUe(u)->m_rnti, u);

        for (size_t g = 0; g < bareGroups.size(); ++g)
        {
            std::vector<UePtrAndBufferReq> grp;
            grp.reserve(bareGroups[g].size());
            for (const auto& uePtr : bareGroups[g])
            {
                grp.push_back(byRnti.at(uePtr->m_rnti));
                rntiToGroup[uePtr->m_rnti] = g;
            }
            groups.push_back(std::move(grp));
        }
    }

    for (size_t g = 0; g < groups.size(); ++g)
    {
        std::string members;
        for (const auto& m : groups[g])
        {
            members += std::to_string(GetUe(m)->m_rnti) + " ";
        }
        NS_LOG_UNCOND("Group " << g << " : " << members);
    }

    // ============================================================
    // PHASE 2: distribute the symbols across the fixed groups.
    // ============================================================
    while (resources > 0)
    {
        if (m_activeDlAi || m_activeUlAi)
        {
            callNotifyFn(ueVector);
        }

        SortUeVector(&ueVector, GetCompareFn);

        long chosenGroup = -1;
        for (auto it = ueVector.begin(); it != ueVector.end(); ++it)
        {
            if (NeedsMoreResources(*it))
            {
                continue;
            }
            chosenGroup = static_cast<long>(rntiToGroup.at(GetUe(*it)->m_rnti));
            break;
        }

        if (chosenGroup < 0)
        {
            NS_LOG_INFO("All groups satisfied; " << resources << " symbols left unused.");
            break;
        }

        auto& group = groups[static_cast<size_t>(chosenGroup)];

        assigned.m_rbg += numOfAssignableRbgs;
        assigned.m_sym += 1;

        for (auto& member : group)
        {
            auto uePtr = GetUe(member);

            auto& assignedRbgs = GetRBGFn(uePtr);
            const auto existingRbgs = assignedRbgs.size();
            assignedRbgs.resize(existingRbgs + numOfAssignableRbgs);
            std::copy(remainingRbgSet.begin(), remainingRbgSet.end(),
                      assignedRbgs.begin() + existingRbgs);

            auto& assignedSymbols = GetSymFn(uePtr);
            const auto existingSymbols = assignedSymbols.size();
            assignedSymbols.resize(existingSymbols + numOfAssignableRbgs);
            std::fill(assignedSymbols.begin() + existingSymbols, assignedSymbols.end(), resources);

            NS_LOG_DEBUG("Assigned " << numOfAssignableRbgs << " " << type
                                     << " RBG to UE " << uePtr->m_rnti
                                     << " (group " << chosenGroup << ", size " << group.size()
                                     << ") label = " << resources);

            SuccessfulAssignmentFn(member, FTResources(numOfAssignableRbgs, 1), assigned);
        }

        for (auto& ue : ueVector)
        {
            if (static_cast<long>(rntiToGroup.at(GetUe(ue)->m_rnti)) != chosenGroup)
            {
                UnSuccessfulAssignmentFn(ue, FTResources(numOfAssignableRbgs, 1), assigned);
            }
        }

        resources -= 1;
    }

    NrMacSchedulerTdmaSdma::BeamSymbolMap ret;
    for (const auto& el : activeUe)
    {
        std::unordered_set<uint32_t> uniqueLabels;
        for (const auto& ue : el.second)
        {
            auto& symVec = GetSymFn(ue.first);
            for (auto lbl : symVec)
            {
                uniqueLabels.insert(static_cast<uint32_t>(lbl));
            }
        }
        ret.insert(std::make_pair(el.first, static_cast<uint32_t>(uniqueLabels.size())));
    }
    return ret;
}

/**
 * @brief Assign the available DL RBG to the UEs
 * @param symAvail Number of available symbols
 * @param activeDl active DL flows and UE
 * @return a map between the beam and the symbols assigned to each one
 *
 * The function will prepare all the needed callbacks to return UE DL parameters
 * (e.g., the DL TBS, the DL RBG) and then will call NrMacSchedulerTdma::AssignRBGTDMA.
 */
NrMacSchedulerTdmaSdma::BeamSymbolMap
NrMacSchedulerTdmaSdma::AssignDLRBG(uint32_t symAvail, const ActiveUeMap& activeDl) const
{
    NS_LOG_FUNCTION(this);
    
    BeforeSchedFn beforeSched = std::bind(&NrMacSchedulerTdmaSdma::BeforeDlSched,
                                          this,
                                          std::placeholders::_1,
                                          std::placeholders::_2);
    AfterSuccessfulAssignmentFn SuccFn = std::bind(&NrMacSchedulerTdmaSdma::AssignedDlResources,
                                                   this,
                                                   std::placeholders::_1,
                                                   std::placeholders::_2,
                                                   std::placeholders::_3);
    AfterUnsuccessfulAssignmentFn UnSuccFn = std::bind(&NrMacSchedulerTdmaSdma::NotAssignedDlResources,
                                                       this,
                                                       std::placeholders::_1,
                                                       std::placeholders::_2,
                                                       std::placeholders::_3);
    GetCompareUeFn compareFn = std::bind(&NrMacSchedulerTdmaSdma::GetUeCompareDlFn, this);

    GetTBSFn GetTbs = &NrMacSchedulerUeInfo::GetDlTBS;
    GetRBGFn GetRBG = &NrMacSchedulerUeInfo::GetDlRBG;
    GetSymFn GetSym = &NrMacSchedulerUeInfo::GetDlSym;

    CallNotifyFn callNotifyFn =
        std::bind(&NrMacSchedulerTdmaSdma::CallNotifyDlFn, this, std::placeholders::_1);

    return AssignRBGTDMA(symAvail,
                         activeDl,
                         "DL",
                         beforeSched,
                         compareFn,
                         GetTbs,
                         GetRBG,
                         GetSym,
                         SuccFn,
                         UnSuccFn,
                         callNotifyFn);
}

/**
 * @brief Assign the available UL RBG to the UEs
 * @param symAvail Number of available symbols
 * @param activeUl active DL flows and UE
 * @return a map between the beam and the symbols assigned to each one
 *
 * The function will prepare all the needed callbacks to return UE UL parameters
 * (e.g., the UL TBS, the UL RBG) and then will call NrMacSchedulerTdma::AssignRBGTDMA.
 */
NrMacSchedulerTdmaSdma::BeamSymbolMap
NrMacSchedulerTdmaSdma::AssignULRBG(uint32_t symAvail, const ActiveUeMap& activeUl) const
{
    NS_LOG_FUNCTION(this);
    BeforeSchedFn beforeSched = std::bind(&NrMacSchedulerTdmaSdma::BeforeUlSched,
                                          this,
                                          std::placeholders::_1,
                                          std::placeholders::_2);
    AfterSuccessfulAssignmentFn SuccFn = std::bind(&NrMacSchedulerTdmaSdma::AssignedUlResources,
                                                   this,
                                                   std::placeholders::_1,
                                                   std::placeholders::_2,
                                                   std::placeholders::_3);
    GetCompareUeFn compareFn = std::bind(&NrMacSchedulerTdmaSdma::GetUeCompareUlFn, this);
    AfterUnsuccessfulAssignmentFn UnSuccFn = std::bind(&NrMacSchedulerTdmaSdma::NotAssignedUlResources,
                                                       this,
                                                       std::placeholders::_1,
                                                       std::placeholders::_2,
                                                       std::placeholders::_3);
    GetTBSFn GetTbs = &NrMacSchedulerUeInfo::GetUlTBS;
    GetRBGFn GetRBG = &NrMacSchedulerUeInfo::GetUlRBG;
    GetSymFn GetSym = &NrMacSchedulerUeInfo::GetUlSym;

    CallNotifyFn callNotifyFn =
        std::bind(&NrMacSchedulerTdmaSdma::CallNotifyUlFn, this, std::placeholders::_1);

    return AssignRBGTDMA(symAvail,
                         activeUl,
                         "UL",
                         beforeSched,
                         compareFn,
                         GetTbs,
                         GetRBG,
                         GetSym,
                         SuccFn,
                         UnSuccFn,
                         callNotifyFn);
}

/**
 * @brief Create a DL DCI starting from spoint and spanning maxSym symbols
 * @param spoint Starting point of the DCI
 * @param ueInfo UE representation
 * @param maxSym Maximum number of symbols for the creation of the DCI
 * @return a pointer to the newly created DCI
 *
 * The method calculates the TBS and the real number of symbols needed, and
 * then call CreateDci().
 */
std::shared_ptr<DciInfoElementTdma>
NrMacSchedulerTdmaSdma::CreateDlDci(PointInFTPlane* spoint,
                                const std::shared_ptr<NrMacSchedulerUeInfo>& ueInfo,
                                [[maybe_unused]] uint32_t maxSym, bool advanceSpoint) const
{
    NS_LOG_FUNCTION(this);
    uint32_t tbs = m_dlAmc->CalculateTbSize(ueInfo->GetDlMcs(),
    // uint32_t tbs = m_dlAmc->CalculateTbSize(12,
                                            ueInfo->m_dlRank,
                                            ueInfo->m_dlRBG.size() * GetNumRbPerRbg());

    // If is less than 10 (3 mac header, 2 rlc header, 5 data), then we can't
    // transmit any new data, so don't create dci.
    if (tbs < 10)
    {
        NS_LOG_DEBUG("While creating DL DCI for UE " << ueInfo->m_rnti << " assigned "
                                                     << ueInfo->m_dlRBG.size()
                                                     << " DL RBG, but TBS < 10");
        ueInfo->m_dlTbSize = 0;
        return nullptr;
    }

    const std::vector<bool> notchedRBGsMask = GetDlNotchedRbgMask();
    int zeroes = std::count(notchedRBGsMask.begin(), notchedRBGsMask.end(), 0);
    uint32_t numOfAssignableRbgs = GetBandwidthInRbg() - zeroes;

    uint8_t numSym = static_cast<uint8_t>(ueInfo->m_dlRBG.size() / numOfAssignableRbgs);

    // std::cout << "NUM SYM = " << +numSym << std::endl;

    auto dci = CreateDci(spoint,
                         ueInfo,
                         tbs,
                         DciInfoElementTdma::DL,
                         ueInfo->GetDlMcs(),
                        //  12,
                         ueInfo->m_dlRank,
                         ueInfo->m_dlPrecMats,
                         std::max(numSym, static_cast<uint8_t>(1)));

    // The starting point must advance.
    if (advanceSpoint)
    {
        spoint->m_rbg = 0;
        spoint->m_sym += numSym;
    }

    return dci;
}

std::shared_ptr<DciInfoElementTdma>
NrMacSchedulerTdmaSdma::CreateDlDci(PointInFTPlane* spoint,
                                const std::shared_ptr<NrMacSchedulerUeInfo>& ueInfo,
                                uint32_t maxSym) const
{
    return CreateDlDci(spoint, ueInfo, maxSym, false);
}

void
NrMacSchedulerTdmaSdma::SortUeVector(std::vector<UePtrAndBufferReq>* ueVector,
                                 const GetCompareUeFn& GetCompareFn) const
{
    std::stable_sort(ueVector->begin(), ueVector->end(), GetCompareFn());
}

/**
 * @brief Create a UL DCI starting from spoint and spanning maxSym symbols
 * @param spoint Starting point of the DCI
 * @param ueInfo UE representation
 * @return a pointer to the newly created DCI
 *
 * The method calculates the TBS and the real number of symbols needed, and
 * then call CreateDci().
 * Allocate the DCI going backward from the starting point (it should be called
 * ending point maybe).
 */
std::shared_ptr<DciInfoElementTdma>
NrMacSchedulerTdmaSdma::CreateUlDci(NrMacSchedulerNs3::PointInFTPlane* spoint,
                                const std::shared_ptr<NrMacSchedulerUeInfo>& ueInfo,
                                uint32_t maxSym) const
{
    NS_LOG_FUNCTION(this);
    uint32_t tbs = m_ulAmc->CalculateTbSize(ueInfo->m_ulMcs,
                                            ueInfo->m_ulRank,
                                            ueInfo->m_ulRBG.size() * GetNumRbPerRbg());

    // If is less than 12, 7 (3 mac header, 2 rlc header, 2 data) + SHORT_BSR (5),
    // then we can't transmit any new data, so don't create dci.
    if (tbs < 12)
    {
        NS_LOG_DEBUG("While creating UL DCI for UE " << ueInfo->m_rnti << " assigned "
                                                     << ueInfo->m_ulRBG.size()
                                                     << " UL RBG, but TBS " << tbs << " < 12");
        return nullptr;
    }

    const std::vector<bool> notchedRBGsMask = GetUlNotchedRbgMask();
    int zeroes = std::count(notchedRBGsMask.begin(), notchedRBGsMask.end(), 0);
    uint32_t numOfAssignableRbgs = GetBandwidthInRbg() - zeroes;

    uint8_t numSym =
        static_cast<uint8_t>(std::max<size_t>(ueInfo->m_ulRBG.size() / numOfAssignableRbgs, 1));
    numSym = std::min(numSym, static_cast<uint8_t>(maxSym));

    NS_ASSERT_MSG(spoint->m_sym >= numSym, "number of symbols are = " << +numSym << ". However, spoint symbols are in " << +spoint->m_sym);

    // The starting point must go backward to accommodate the needed sym
    spoint->m_sym -= numSym;

    auto dci = CreateDci(spoint,
                         ueInfo,
                         tbs,
                         DciInfoElementTdma::UL,
                         ueInfo->m_ulMcs,
                         ueInfo->m_ulRank,
                         ueInfo->m_ulPrecMats,
                         numSym);

    // Reset the RBG (we are TDMA)
    spoint->m_rbg = 0;

    return dci;
}

uint8_t
NrMacSchedulerTdmaSdma::GetTpc() const
{
    NS_LOG_FUNCTION(this);
    return 1; // 1 is mapped to 0 for Accumulated mode, and to -1 in Absolute mode TS38.213 Table
              // Table 7.1.1-1
}


/**
 * @brief Create a DCI with the parameters specified as input
 * @param spoint starting point
 * @param ueInfo ue representation
 * @param tbs Transport Block Size
 * @param fmt Format of the DCI (UL or DL)
 * @param mcs MCS
 * @param numSym Number of symbols
 * @return a pointer to the newly created DCI
 *
 * Creates a TDMA DCI (a DCI with all the resource block assigned for the
 * specified number of symbols)
 */
std::shared_ptr<DciInfoElementTdma>
NrMacSchedulerTdmaSdma::CreateDci(NrMacSchedulerNs3::PointInFTPlane* spoint,
                              const std::shared_ptr<NrMacSchedulerUeInfo>& ueInfo,
                              uint32_t tbs,
                              DciInfoElementTdma::DciFormat fmt,
                              uint32_t mcs,
                              uint8_t rank,
                              Ptr<const ComplexMatrixArray> precMats,
                              uint8_t numSym) const
{
    NS_LOG_FUNCTION(this);
    NS_ASSERT(tbs > 0);
    NS_ASSERT(numSym > 0);

    std::shared_ptr<DciInfoElementTdma> dci =
        std::make_shared<DciInfoElementTdma>(ueInfo->m_rnti,
                                             fmt,
                                             spoint->m_sym,
                                             numSym,
                                             mcs,
                                             rank,
                                             precMats,
                                             tbs,
                                             1,
                                             0,
                                             DciInfoElementTdma::DATA,
                                             GetBwpId(),
                                             GetTpc());

    std::vector<bool> rbgAssigned =
        fmt == DciInfoElementTdma::DL ? GetDlNotchedRbgMask() : GetUlNotchedRbgMask();

    if (rbgAssigned.empty())
    {
        rbgAssigned = std::vector<bool>(GetBandwidthInRbg(), true);
    }

    NS_ASSERT(rbgAssigned.size() == GetBandwidthInRbg());

    dci->m_rbgBitmask = std::move(rbgAssigned);

    std::ostringstream oss;
    for (auto x : dci->m_rbgBitmask)
    {
        oss << std::to_string(x) << " ";
    }

    NS_LOG_UNCOND("TIME: " << Simulator::Now().GetMicroSeconds() << " UE " << ueInfo->m_rnti << " assigned RBG from " << spoint->m_rbg << " with mask "
                      << oss.str() << " for " << static_cast<uint32_t>(numSym) << " SYM " << "in symbol " << +spoint->m_sym);

    NS_ASSERT(std::count(dci->m_rbgBitmask.begin(), dci->m_rbgBitmask.end(), 0) !=
              GetBandwidthInRbg());

    return dci;
}


std::vector<DciInfoElementTdma>
NrMacSchedulerTdmaSdma::DoReshapeAllocation(
    const std::vector<DciInfoElementTdma>& dcis,
    uint8_t& startingSymbol,
    uint8_t& numSymbols,
    std::vector<bool>& bitmask,
    const bool isDl,
    const std::unordered_map<uint16_t, std::shared_ptr<NrMacSchedulerUeInfo>>& ueMap)
{
    // Declare lambda function to compute MCS based on sub-band information
    auto computeMcs = [isDl](const std::shared_ptr<NrMacSchedulerUeInfo>& ueInfo,
                             std::vector<uint16_t>& rbgVector) -> double {
        if (isDl)
        {
            double currentMcs = ueInfo->m_dlMcs; // Wideband MCS
            if (!rbgVector.empty())
            {
                const auto sum = std::transform_reduce(
                    rbgVector.begin(),
                    rbgVector.end(),
                    0.0,
                    std::plus<>(),
                    [ueInfo](auto a) {
                        return ueInfo->m_dlSbMcsInfo.at(ueInfo->m_rbgToSb.at(a)).mcs;
                    });
                currentMcs = sum / rbgVector.size();
            }
            return currentMcs;
        }
        else
        {
            return ueInfo->m_ulMcs;
        }
    };
    // clang-format off
    /**
     * TDMA DCI consolidation/defragmentation follows these steps
     * 1. Pick a DCI
     * 2. Compute number of resources required by DCI
     * 3. Give all available RBGs to UE
     * 4. Sort RBGs based on best sub-band
     * 5. Check if we have any chance of meeting the number of resources at the MCS specified at the DCI
     * 5.1 If not, we try to remove the lowest RBG (go back to 4).
     * 5.2 If yes, we found our allocation, continue.
     * 6. If this is not the last DCI and there are remaining RBGs, go back to 1. Else, continue.
     * 7. If all DCIs were allocated, and we still have RBGs available, try to reduce number of
     * symbols used, by spreading DCIs in remaining RBGs, to free up symbols to other beams.
     */
    // clang-format on
    uint8_t availableSymbols = numSymbols;
    std::vector<DciInfoElementTdma> reshapedDcis{};

    // Step 1, pick a dci
    for (auto dci : dcis)
    {
        auto& ueInfo = ueMap.at(dci.m_rnti);

        // Step 2, compute the number of resources needed
        const std::size_t numResources =
            dci.m_numSym * std::count(dci.m_rbgBitmask.begin(), dci.m_rbgBitmask.end(), true);

        // Step 3, allocate all RBGs to UE
        std::vector<uint16_t> allocatedRbgs;
        for (std::size_t i = 0; i < bitmask.size(); i++)
        {
            if (bitmask.at(i))
            {
                allocatedRbgs.push_back(i);
            }
        }

        // We want to find the set of RBGs that return the maximum MCS
        if (isDl && !ueInfo->m_rbgToSb.empty())
        {
            auto prevMcs = computeMcs(ueInfo, allocatedRbgs);

            // Step 4, sort RBGs based on sub-band MCS (from highest to lowest)
            std::stable_sort(allocatedRbgs.begin(),
                             allocatedRbgs.end(),
                             [&](uint16_t a, uint16_t b) {
                                 return ueInfo->m_dlSbMcsInfo.at(ueInfo->m_rbgToSb.at(a)).mcs >
                                        ueInfo->m_dlSbMcsInfo.at(ueInfo->m_rbgToSb.at(b)).mcs;
                             });
            // While we have remaining RBGs and the DCI number of resources fit into the remaining
            // resources, we try to remove bad RBGs to increase overall MCS
            while ((!allocatedRbgs.empty()) &&
                   (((allocatedRbgs.size() - 1) * availableSymbols) >= numResources))
            {
                if (ueInfo->m_dlSbMcsInfo.at(ueInfo->m_rbgToSb.at(allocatedRbgs.back())).mcs >=
                    std::min(dci.m_mcs, ueInfo->m_dlMcs))
                {
                    // There will be no MCS improvement in removing additional RBGs
                    break;
                }

                // If sub-band MCS is lower than wideband, we take this RBG out
                auto currRbg = allocatedRbgs.back();
                allocatedRbgs.pop_back();
                auto currMcs = computeMcs(ueInfo, allocatedRbgs);

                // Things will only get worse if we continue removing
                if (currMcs <= prevMcs)
                {
                    allocatedRbgs.push_back(currRbg);
                    break;
                }
            }
        }

        // Compute the number of required symbols
        uint8_t minSymbols = ceil((double)numResources / allocatedRbgs.size());

        // Remove RBGs or increase number of symbols until we match the number of resources
        uint32_t currResources = minSymbols * allocatedRbgs.size();
        while ((numResources != currResources) && (!allocatedRbgs.empty()))
        {
            if (currResources < numResources)
            {
                minSymbols += 1;
            }
            else if (currResources > numResources)
            {
                allocatedRbgs.pop_back();
            }
            currResources = minSymbols * allocatedRbgs.size();
        }

        if (minSymbols <= availableSymbols && currResources == numResources)
        {
            std::vector<bool> allocatedBitmask(bitmask.size(), false);
            for (auto rbg : allocatedRbgs)
            {
                allocatedBitmask.at(rbg) = true;
            }
            // Update DCI after reshaping
            reshapedDcis.emplace_back(startingSymbol, minSymbols, allocatedBitmask, dci);
            availableSymbols -= minSymbols;
            startingSymbol += minSymbols;
            numSymbols -= minSymbols;
        }
    }
    return reshapedDcis;
}

// =============================================================================
//  (A) INTEGRATION POINT: real DL channel matrix for one UE [Nrx x Ntx x Pages]
// =============================================================================
Ptr<const ComplexMatrixArray>
NrMacSchedulerTdmaSdma::GetDlChannelMatrix(
    const std::shared_ptr<NrMacSchedulerUeInfo>& ueInfo) const
{
    NS_LOG_FUNCTION(this << ueInfo->m_rnti);
    // TODO(you): return ueInfo's cached DL channel, e.g. ueInfo->m_dlChannelMatrix,
    // using the same source ComputeGroupZfPrecoding() reads from.
    
    return ueInfo->m_dlChanMat;
}

// =============================================================================
//  (E) Effective per-layer spectral efficiency -> MCS index.
//      Replace the table with your AMC mapping for fidelity.
// =============================================================================
uint8_t
SeToMcs(double sePerLayer)
{
    // --- REPLACE WITH YOUR AMC, e.g. -----------------------------------------
    //   return m_dlAmc->GetMcsFromSpectralEfficiency(sePerLayer);   // (verify name)
    // -------------------------------------------------------------------------
    // Fallback: NR 256-QAM table (TS 38.214, table 5.1.3.1-2) SE thresholds.
    static const double kSe[] = {
        0.2344, 0.3770, 0.6016, 0.8770, 1.1758, 1.4766, 1.6953, 1.9141, 2.1602,
        2.4063, 2.5703, 2.7305, 3.0293, 3.3223, 3.6094, 3.9023, 4.2129, 4.5234,
        4.8164, 5.1152, 5.3320, 5.5547, 5.8906, 6.2266, 6.5703, 6.9141, 7.1602, 7.4063};
    uint8_t mcs = 0;
    NS_LOG_UNCOND("SePerLayer: " << sePerLayer);
    for (uint8_t i = 0; i < sizeof(kSe) / sizeof(kSe[0]); ++i)
        if (sePerLayer >= kSe[i]) mcs = i;
    return mcs;
}

// =============================================================================
//  Commit precoders + rank, then recompute MCS from the effective channel.
//  Interference is summed over EVERY transmitted column of EVERY UE, so:
//    - ZF        -> interference ~ 0  (SINR = SNR, with ZF noise enhancement
//                                      already baked into the column norms)
//    - aggregate -> real inter-user interference, lowering the MCS
// =============================================================================
void
FinalizeAllocation(std::vector<UeAlloc>& allocs,
                                           uint32_t Ntx,
                                           size_t pages,
                                           double noiseLin,
                                           double sinrBackOffDb,
                                           int forceDlMcs)
{
    // 1) write per-UE precoder ComplexMatrixArray + rank
    for (auto& a : allocs)
    {
        if (a.rank == 0)
        {
            a.ue->m_dlRank = 0;
            a.ue->m_dlPrecMats = nullptr;
            continue;
        }
        Ptr<ComplexMatrixArray> P = Create<ComplexMatrixArray>(Ntx, a.rank, pages); // (B)
        for (size_t p = 0; p < pages; ++p)
            for (uint32_t r = 0; r < Ntx; ++r)
                for (uint32_t c = 0; c < a.rank; ++c)
                    (*P)(r, c, p) = a.Wp[p](r, c);                                   // (B)
        a.ue->m_dlRank = static_cast<uint8_t>(a.rank);
        a.ue->m_dlPrecMats = P;
    }

    // 2) per-UE effective SINR -> per-layer SE -> MCS
    const double N0 = (noiseLin > 0.0) ? noiseLin : 1.0;          // (D)
    const double sinrScale = std::pow(10.0, -sinrBackOffDb / 10.0);

    for (auto& a : allocs)
    {
        if (a.rank == 0)
            continue;

        double seAccum = 0.0;     // sum of log2(1+SINR) over (streams, pages)
        uint32_t terms = 0;

        for (size_t p = 0; p < pages; ++p)
        {
            const Eigen::MatrixXcd Hk = PageToEigen(a.H, p);   // Nrx x Ntx

            for (uint32_t i = 0; i < a.rank; ++i)
            {
                const Eigen::VectorXcd u = a.Up[p].col(i);      // combiner
                const Eigen::RowVectorXcd uH_Hk = u.adjoint() * Hk;  // 1 x Ntx

                const std::complex<double> sig = uH_Hk * a.Wp[p].col(i);
                const double signal = std::norm(sig);

                double interf = 0.0;
                for (const auto& b : allocs)               // every other column
                {
                    if (b.rank == 0)
                        continue;
                    for (uint32_t j = 0; j < b.rank; ++j)
                    {
                        if (&b == &a && j == i)
                            continue;
                        const std::complex<double> x = uH_Hk * b.Wp[p].col(j);
                        interf += std::norm(x);
                    }
                }

                const double sinr = signal / (interf + N0) * sinrScale;
                seAccum += std::log2(1.0 + sinr);
                ++terms;
            }
        }

        const double sePerLayer = (terms > 0) ? seAccum / terms : 0.0;
        uint8_t mcs = SeToMcs(sePerLayer);                                     // (E)
        if (forceDlMcs >= 0) mcs = static_cast<uint8_t>(std::min(forceDlMcs, 27));

        a.ue->m_dlMcs = mcs;                                                         // (F)

        // NS_LOG_DEBUG("UE " << a.ue->m_rnti << " rank=" << a.rank
        //                    << " SE/layer=" << sePerLayer << " -> MCS " << +mcs);
    }
}

// =============================================================================
//  ZF / eigen-ZF baseline (AggregateSvd = false)
// =============================================================================
void
NrMacSchedulerTdmaSdma::ComputeFullCsiSvdPrecoding(const std::vector<UePtr>& group,
                                                   uint32_t maxStreams) const
{
    NS_LOG_FUNCTION(this << group.size() << maxStreams);

    struct Csi
    {
        UePtr ue;
        Ptr<const ComplexMatrixArray> H;
        std::vector<Eigen::MatrixXcd> V, U;   // per page
        std::vector<Eigen::VectorXd> sv;
        uint32_t rank = 0;
    };

    std::vector<Csi> csi;
    size_t pages = 0;
    uint32_t Ntx = 0;
    for (const auto& ue : group)
    {
        if (!ue) continue;
        auto H = GetDlChannelMatrix(ue);
        if (!H || H->GetNumPages() == 0 || H->GetNumCols() == 0) continue;
        Csi e; e.ue = ue; e.H = H;
        pages = H->GetNumPages(); Ntx = H->GetNumCols();
        csi.push_back(std::move(e));
    }
    if (csi.empty() || Ntx == 0 || pages == 0) return;

    for (auto& e : csi)
    {
        e.V.resize(pages); e.U.resize(pages); e.sv.resize(pages);
        for (size_t p = 0; p < pages; ++p)
        {
            Eigen::MatrixXcd Hk = PageToEigen(e.H, p);
            Eigen::JacobiSVD<Eigen::MatrixXcd> svd(Hk, Eigen::ComputeThinU | Eigen::ComputeThinV);
            e.U[p] = svd.matrixU();
            e.V[p] = svd.matrixV();
            e.sv[p] = svd.singularValues();
        }
    }

    // choose streams (rank) honoring the cap
    struct Cand { size_t k; uint32_t s; double score; };
    std::vector<Cand> cands;
    for (size_t k = 0; k < csi.size(); ++k)
        for (uint32_t s = 0; s < static_cast<uint32_t>(csi[k].sv[0].size()); ++s)
        {
            double acc = 0.0;
            for (size_t p = 0; p < pages; ++p) acc += csi[k].sv[p](s);
            cands.push_back({k, s, acc / pages});
        }
    if (cands.empty()) return;
    std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b){ return a.score > b.score; });

    uint32_t cap = (maxStreams == 0) ? Ntx : std::min(maxStreams, Ntx);
    cap = std::min<uint32_t>(cap, static_cast<uint32_t>(cands.size()));
    

    // which (k,s) streams are kept, in priority order
    // std::vector<std::pair<size_t,uint32_t>> kept;
    // std::vector<uint32_t> keepCount(csi.size(), 0);
    // for (uint32_t i = 0; i < cap; ++i) { kept.push_back({cands[i].k, cands[i].s}); keepCount[cands[i].k]++; }
    // for (size_t k = 0; k < csi.size(); ++k) csi[k].rank = keepCount[k];

    const uint32_t maxRankPerUe = std::max<uint32_t>(m_maxRankPerUe, 1);   // NEW
    std::vector<std::pair<size_t,uint32_t>> kept;
    std::vector<uint32_t> keepCount(csi.size(), 0);
    for (size_t i = 0; i < cands.size() && kept.size() < cap; ++i)
    {
        const size_t k = cands[i].k;
        if (keepCount[k] >= maxRankPerUe) continue;        // <-- per-UE rank cap
        kept.push_back({cands[i].k, cands[i].s});
        keepCount[k]++;
    }
    for (size_t k = 0; k < csi.size(); ++k) csi[k].rank = keepCount[k];

    uint32_t S = 0; for (auto& e : csi) S += e.rank;
    if (S == 0) return;

    // build per-UE allocation
    std::vector<UeAlloc> allocs(csi.size());
    for (size_t k = 0; k < csi.size(); ++k)
    {
        allocs[k].ue = csi[k].ue; allocs[k].H = csi[k].H; allocs[k].rank = csi[k].rank;
        allocs[k].Wp.assign(pages, Eigen::MatrixXcd(Ntx, std::max<uint32_t>(csi[k].rank,1)));
        allocs[k].Up.assign(pages, Eigen::MatrixXcd());
    }

    for (size_t p = 0; p < pages; ++p)
    {
        // stack kept unit directions -> G (S x Ntx)
        Eigen::MatrixXcd G(S, Ntx);
        for (uint32_t r = 0; r < S; ++r)
            G.row(r) = csi[kept[r].first].V[p].col(kept[r].second).adjoint();

        Eigen::MatrixXcd W = G.completeOrthogonalDecomposition().pseudoInverse(); // Ntx x S, G W = I

        // unit-norm columns, equal power split (ZF separation preserved)
        for (uint32_t c = 0; c < S; ++c) { double n = W.col(c).norm(); if (n > 0) W.col(c) /= n; }
        W /= std::sqrt(static_cast<double>(S));

        // split columns + matching combiners to the owning UE
        std::vector<uint32_t> localCol(csi.size(), 0);
        for (uint32_t r = 0; r < S; ++r)
        {
            const size_t k = kept[r].first;
            const uint32_t s = kept[r].second;
            const uint32_t lc = localCol[k]++;
            allocs[k].Wp[p].col(lc) = W.col(r);
            if (allocs[k].Up[p].size() == 0)
                allocs[k].Up[p].resize(csi[k].U[p].rows(), csi[k].rank);
            allocs[k].Up[p].col(lc) = csi[k].U[p].col(s);   // left singular vec = combiner
        }
    }

    FinalizeAllocation(allocs, Ntx, pages, m_dlNoisePowerLin, m_sinrBackoffDb, m_forceDlMcs);
}

// =============================================================================
//  Aggregate-SVD baseline (AggregateSvd = true)
//  Stack ALL channels, SVD, use right singular vectors as a joint precoder,
//  assign each eigen-stream to the UE that receives it strongest.
// =============================================================================
void
NrMacSchedulerTdmaSdma::ComputeAggregateSvdPrecoding(const std::vector<UePtr>& group,
                                                     uint32_t maxStreams) const
{
    NS_LOG_FUNCTION(this << group.size() << maxStreams);

    std::vector<UePtr> ues;
    std::vector<Ptr<const ComplexMatrixArray>> H;
    size_t pages = 0; uint32_t Ntx = 0;
    for (const auto& ue : group)
    {
        if (!ue) continue;
        auto h = GetDlChannelMatrix(ue);
        if (!h || h->GetNumPages() == 0 || h->GetNumCols() == 0) continue;
        ues.push_back(ue); H.push_back(h);
        pages = h->GetNumPages(); Ntx = h->GetNumCols();
    }
    if (ues.empty() || Ntx == 0 || pages == 0) return;

    // per-page aggregate SVD; cache V and per-UE channel pages
    std::vector<Eigen::MatrixXcd> Vp(pages);                       // Ntx x nSV
    std::vector<std::vector<Eigen::MatrixXcd>> Hk(ues.size());     // [ue][page]
    for (size_t i = 0; i < ues.size(); ++i) Hk[i].resize(pages);

    uint32_t nSV = 0;
    for (size_t p = 0; p < pages; ++p)
    {
        uint32_t R = 0;
        for (size_t i = 0; i < ues.size(); ++i) { Hk[i][p] = PageToEigen(H[i], p); R += Hk[i][p].rows(); }
        Eigen::MatrixXcd Hagg(R, Ntx);
        uint32_t row = 0;
        for (size_t i = 0; i < ues.size(); ++i) { Hagg.middleRows(row, Hk[i][p].rows()) = Hk[i][p]; row += Hk[i][p].rows(); }
        Eigen::JacobiSVD<Eigen::MatrixXcd> svd(Hagg, Eigen::ComputeThinV);
        Vp[p] = svd.matrixV();
        nSV = std::max<uint32_t>(nSV, static_cast<uint32_t>(svd.matrixV().cols()));
    }

    uint32_t cap = (maxStreams == 0) ? std::min(nSV, Ntx) : std::min({maxStreams, nSV, Ntx});
    if (cap == 0) return;

    // assign each of the top-`cap` eigen-streams to the UE that receives it strongest
    // (decided wideband: average ||H_k v_c||^2 over pages)
    std::vector<int> owner(cap, -1);
    for (uint32_t c = 0; c < cap; ++c)
    {
        double best = -1.0; int bestUe = -1;
        for (size_t i = 0; i < ues.size(); ++i)
        {
            double e = 0.0;
            for (size_t p = 0; p < pages; ++p) { if (c < Vp[p].cols()) e += (Hk[i][p] * Vp[p].col(c)).squaredNorm(); }
            if (e > best) { best = e; bestUe = static_cast<int>(i); }
        }
        owner[c] = bestUe;
    }

    std::vector<uint32_t> rank(ues.size(), 0);
    for (uint32_t c = 0; c < cap; ++c) rank[owner[c]]++;

    std::vector<UeAlloc> allocs(ues.size());
    for (size_t i = 0; i < ues.size(); ++i)
    {
        allocs[i].ue = ues[i]; allocs[i].H = H[i]; allocs[i].rank = rank[i];
        allocs[i].Wp.assign(pages, Eigen::MatrixXcd(Ntx, std::max<uint32_t>(rank[i],1)));
        allocs[i].Up.assign(pages, Eigen::MatrixXcd(Hk[i][0].rows(), std::max<uint32_t>(rank[i],1)));
    }

    const double invSqrtS = 1.0 / std::sqrt(static_cast<double>(cap)); // equal power split
    std::vector<uint32_t> localCol(ues.size(), 0);
    // assign columns in stream order so localCol mapping is consistent across pages
    std::vector<uint32_t> colOfStream(cap);
    { std::vector<uint32_t> lc(ues.size(),0); for (uint32_t c=0;c<cap;++c){ colOfStream[c]=lc[owner[c]]++; } }

    for (size_t p = 0; p < pages; ++p)
        for (uint32_t c = 0; c < cap; ++c)
        {
            const int k = owner[c];
            const uint32_t lc = colOfStream[c];
            // Eigen::VectorXcd v = (c < Vp[p].cols()) ? Vp[p].col(c) : Eigen::VectorXcd::Zero(Ntx);
            Eigen::VectorXcd v;
            if (c < static_cast<uint32_t>(Vp[p].cols())) v = Vp[p].col(c);
            else                                         v = Eigen::VectorXcd::Zero(Ntx);
            allocs[k].Wp[p].col(lc) = v * invSqrtS;            // joint precoder column
            Eigen::VectorXcd rx = Hk[k][p] * v;                 // matched-filter combiner
            double n = rx.norm(); if (n > 0) rx /= n;
            allocs[k].Up[p].col(lc) = rx;
        }

    FinalizeAllocation(allocs, Ntx, pages, m_dlNoisePowerLin, m_sinrBackoffDb, m_forceDlMcs);
}

}