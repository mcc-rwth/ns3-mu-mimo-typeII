#include "nr-pm-search-type-two.h"

#include "ns3/boolean.h"
#include "ns3/double.h"
#include "ns3/enum.h"
#include "ns3/simulator.h"

#include "fstream"
#include "complex"
#include "numeric"

#include <iostream>
#include <typeinfo>
#include <sys/stat.h>

#include <vector>
#include <limits>
#include <cstdint>
#include <cmath>
#include <algorithm>

// Eigen library includes
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <unsupported/Eigen/KroneckerProduct>

namespace ns3
{

NS_LOG_COMPONENT_DEFINE("NrPmSearchTypeTwo");
NS_OBJECT_ENSURE_REGISTERED(NrPmSearchTypeTwo);

static bool
FileExists(const std::string& name)
{
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

static void
DumpComplexMatrixArrayToCsv(const ComplexMatrixArray& H,
                            double t,
                            const std::string& filename)
{
    const bool exists = FileExists(filename);

    std::ofstream file(filename, std::ios::out | std::ios::app);
    if (!file.is_open())
    {
        NS_LOG_ERROR("Cannot open file " << filename);
        return;
    }

    if(!exists)
    {
        file << "t,rb,row,col,real,imag\n";
    }

    // header
    // file << "t,rb,row,col,real,imag\n";

    const size_t nRows  = H.GetNumRows();
    const size_t nCols  = H.GetNumCols();
    const size_t nPages = H.GetNumPages(); // RBs

    for (size_t rb = 0; rb < nPages; ++rb)
    {
        for (size_t r = 0; r < nRows; ++r)
        {
            for (size_t c = 0; c < nCols; ++c)
            {
                std::complex<double> v = H(r, c, rb);

                file << t << ","
                     << rb << ","
                     << r  << ","
                     << c  << ","
                     << v.real() << ","
                     << v.imag() << "\n";
            }
        }
    }

    file.close();
}

struct BeamSelResult
    {
        std::vector<uint32_t> m1;
        std::vector<uint32_t> m2;
        std::vector<Eigen::VectorXcd> v;
    };

static inline Eigen::VectorXcd
MakeSteerVec(uint32_t N, uint32_t O, uint32_t idx)
{
    const double denom = static_cast<double>(N) * static_cast<double>(O);
    Eigen::VectorXcd a(N);
    for (uint32_t n = 0; n < N; ++n)
    {
        const double phase = 2 * M_PI * static_cast<double>(idx) * static_cast<double>(n) / denom;
        a(n) = std::exp(std::complex<double>(0.0, phase));
    }

    return a;
}

static inline double
BeamMetric(const Eigen::MatrixXcd& R,
            const Eigen::MatrixXcd& v_l,
            const Eigen::MatrixXcd& u_m)
{
    Eigen::VectorXcd v = Eigen::kroneckerProduct(v_l, u_m);
    const std::complex<double> s = v.adjoint() * (R * v);

    return std::real(s);
}

static inline BeamSelResult
SelectBeamsTypeTwo(const Eigen::MatrixXcd& R,
                    uint32_t N1, uint32_t N2,
                    uint32_t O1, uint32_t O2,
                    uint32_t L)
{
    NS_ASSERT_MSG(R.rows() == R.cols(), "R must be square");

    const uint32_t Nant = static_cast<uint32_t>(R.rows());
    NS_ASSERT_MSG(N1 * N2 == Nant, "Require N1 x N2 == size(R) for kron(v_l, u_m)");

    BeamSelResult out;
    out.m1.reserve(L);
    out.m2.reserve(L);
    out.v.reserve(L);

    // Make a set out of N1 and N2
    std::vector<uint32_t> setN1(N1), setN2(N2);
    for (uint32_t i = 0; i < N1; ++i) setN1[i] = i;
    for (uint32_t i = 0; i < N2; ++i) setN2[i] = i;

    uint32_t q1 = 0, q2 = 0;

    for (uint32_t i = 0; i < L; ++i)
    {
        if (i == 0)
        {
            double best = - std::numeric_limits<double>::infinity();
            uint32_t best_m1 = 0, best_m2 = 0;

            for (uint32_t l= 0; l < N1 * O1; ++l)
            {
                Eigen::VectorXcd v_l = MakeSteerVec(N1, O1, l);
                for (uint32_t m = 0; m < N2 * O2; ++m)
                {
                    Eigen::VectorXcd u_m = MakeSteerVec(N2, O2, m);
                    const double metric = BeamMetric(R, v_l, u_m);
                    if (metric > best)
                    {
                        best = metric;
                        best_m1 = l;
                        best_m2 = m;
                    }
                }
            }

            q1 = best_m1 % O1;
            q2 = best_m2 % O2;

            const uint32_t n1 = best_m1 / O1;
            const uint32_t n2 = best_m2 / O2;

            setN1.erase(std::remove(setN1.begin(), setN1.end(), n1), setN1.end());
            setN2.erase(std::remove(setN2.begin(), setN2.end(), n2), setN2.end());

            // Store the results
            out.m1.push_back(best_m1);
            out.m2.push_back(best_m2);
            out.v.push_back(Eigen::kroneckerProduct(MakeSteerVec(N1, O1, best_m1), MakeSteerVec(N2, O2, best_m2)));
        }
        else
        {
            if (setN1.empty() || setN2.empty())
            {
                break;
            }

            double best = - std::numeric_limits<double>::infinity();
            uint32_t best_n1 = setN1[0];
            uint32_t best_n2 = setN2[0];

            for (uint32_t cand_n1 : setN1)
            {
                const uint32_t l = cand_n1 * O1 + q1;
                Eigen::VectorXcd v_l = MakeSteerVec(N1, O1, l);

                for (uint32_t cand_n2 : setN2)
                {
                    const uint32_t m = cand_n2 * O2 + q2;
                    Eigen::VectorXcd u_m = MakeSteerVec(N2, O2, m);

                    const double metric = BeamMetric(R, v_l, u_m);
                    if (metric > best)
                    {
                        best = metric;
                        best_n1 = cand_n1;
                        best_n2 = cand_n2;
                    }
                }
            }

            const uint32_t sel_m1 = best_n1 * O1 + q1;
            const uint32_t sel_m2 = best_n2 * O2 + q2;

            setN1.erase(std::remove(setN1.begin(), setN1.end(), best_n1), setN1.end());
            setN2.erase(std::remove(setN2.begin(), setN2.end(), best_n2), setN2.end());

            // Store the result
            out.m1.push_back(sel_m1);
            out.m2.push_back(sel_m2);
            out.v.push_back(Eigen::kroneckerProduct(MakeSteerVec(N1, O1, sel_m1), MakeSteerVec(N2, O2, sel_m2)));
        }
    }

    return out;
}


static inline ComplexMatrixArray
RightMultiply(const NrIntfNormChanMat& H,
              const Eigen::MatrixXcd& W1)
{
    const size_t nRows  = H.GetNumRows();   // e.g., nRxPorts
    const size_t nCols  = H.GetNumCols();   // e.g., 32 (=2*Nant)
    const size_t nPages = H.GetNumPages();  // e.g., 7

    NS_ASSERT_MSG(W1.rows() == static_cast<int>(nCols),
                  "RightMultiply: W1 rows must equal H cols");

    const size_t outCols = static_cast<size_t>(W1.cols());

    ComplexMatrixArray HA(nRows, outCols, nPages);

    for (size_t p = 0; p < nPages; ++p)
    {
        // Map H(page) as Eigen view
        Eigen::Map<const Eigen::MatrixXcd> Hpage(
            H.GetPagePtr(p),
            nRows,
            nCols);

        // Map output page memory
        Eigen::Map<Eigen::MatrixXcd> HApage(
            HA.GetPagePtr(p),
            nRows,
            outCols);

        // Multiply
        HApage.noalias() = Hpage * W1;
    }

    return HA;
}


struct W2Result
{
    Eigen::MatrixXcd w2;   // (2L x SB)
    Eigen::MatrixXcd w22;  // (2L x SB)
    Eigen::MatrixXcd W2;   // normalized (2L x SB)
    Eigen::MatrixXcd W22;  // normalized (2L x SB)
    size_t i_pr = 0;
    size_t i_pr2 = 0;
};

static inline W2Result
ComputeW2andNormalize(const ComplexMatrixArray& HA)
{
    const size_t Nr   = HA.GetNumRows();
    const size_t twoL = HA.GetNumCols();
    const size_t SB   = HA.GetNumPages();

    W2Result out;
    out.w2  = Eigen::MatrixXcd::Zero(twoL, SB);
    out.w22 = Eigen::MatrixXcd::Zero(twoL, SB);

    // 1) Per-subband eigvecs of X = HA^H HA
    for (size_t t = 0; t < SB; ++t)
    {
        Eigen::Map<const Eigen::MatrixXcd> HApage(HA.GetPagePtr(t), Nr, twoL);

        Eigen::MatrixXcd X = HApage.adjoint() * HApage;

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(X);
        NS_ASSERT_MSG(es.info() == Eigen::Success, "Eigen decomposition failed");

        const auto& V = es.eigenvectors(); // columns sorted by ascending eigenvalue

        out.w2.col(t)  = V.col(twoL - 1);  // largest
        out.w22.col(t) = V.col(twoL - 2);  // second largest
    }

    // 2) i_pr = argmax_row sum_t |w2(row,t)|^2
    //    i_pr2 = argmax_row sum_t |w22(row,t)|^2
    auto rowEnergyArgmax = [&](const Eigen::MatrixXcd& M) -> size_t
    {
        double best = -1.0;
        size_t bestIdx = 0;
        for (size_t r = 0; r < static_cast<size_t>(M.rows()); ++r)
        {
            // sum |.|^2 over columns
            double e = M.row(r).cwiseAbs2().sum();
            if (e > best)
            {
                best = e;
                bestIdx = r;
            }
        }
        return bestIdx;
    };

    out.i_pr  = rowEnergyArgmax(out.w2);
    out.i_pr2 = rowEnergyArgmax(out.w22);

    // 3) W2 = w2 / w2[i_pr, :]   (elementwise division by the reference row, per column)
    //    W22 = w22 / w22[i_pr2, :]
    out.W2  = out.w2;
    out.W22 = out.w22;

    const double eps = 1e-12; // avoid divide-by-zero blowups

    for (size_t t = 0; t < SB; ++t)
    {
        std::complex<double> denom1 = out.w2(out.i_pr, static_cast<int>(t));
        std::complex<double> denom2 = out.w22(out.i_pr2, static_cast<int>(t));

        if (std::abs(denom1) < eps)
        {
            // If denom is ~0, you must decide policy. Here: leave column unnormalized.
            // Alternative: pick next-best i_pr for that t, or set column to 0.
            continue;
        }
        if (std::abs(denom2) < eps)
        {
            continue;
        }

        out.W2.col(t)  /= denom1;
        out.W22.col(t) /= denom2;
    }

    return out;
}

static const std::vector<double> LEVELS = {
    0.0,
    std::sqrt(1.0/64.0),
    std::sqrt(1.0/32.0),
    std::sqrt(1.0/16.0),
    std::sqrt(1.0/8.0),
    std::sqrt(1.0/4.0),
    std::sqrt(1.0/2.0),
    1.0
};

static std::vector<double> BuildThresholds()
{
    std::vector<double> thr(LEVELS.size() - 1);
    for (size_t i = 0; i < thr.size(); ++i)
    {
        thr[i] = 0.5 * (LEVELS[i] + LEVELS[i+1]);
    }
    return thr;
}

struct QuantResult
{
    Eigen::VectorXi k;      // quantization index
    Eigen::VectorXd p;      // reconstruction level
};

static QuantResult
Q_WB(const Eigen::VectorXd& x)
{
    const auto THR = BuildThresholds();
    const int N = x.size();

    QuantResult out;
    out.k = Eigen::VectorXi::Zero(N);
    out.p = Eigen::VectorXd::Zero(N);

    for (int i = 0; i < N; ++i)
    {
        double xi = std::clamp(x(i), 0.0, 1.0);

        int bin = 0;
        while (bin < static_cast<int>(THR.size()) && xi > THR[bin])
        {
            bin++;
        }

        out.k(i) = bin;            // 0..7
        out.p(i) = LEVELS[bin];
    }

    return out;
}

Eigen::MatrixXcd ComputeD(const Eigen::MatrixXcd& W2,
                           const Eigen::VectorXd& p)
{
    const int rows = W2.rows();
    const int SB   = W2.cols();

    Eigen::MatrixXcd D(rows, SB);

    const double eta = 1e-6;
    double maxp = p.maxCoeff();
    double eps = (maxp > 0.0) ? eta * maxp : 1e-12;

    for (int i = 0; i < rows; ++i)
    {
        double denom = std::max(p(i), eps);
        D.row(i) = W2.row(i) / denom;
    }

    return D;
}

struct QSBResult
{
    Eigen::MatrixXi k2;      // 0 or 1
    Eigen::MatrixXd p2;      // reconstructed magnitude
};

QSBResult Q_SB(const Eigen::MatrixXd& absD)
{
    const double lvl0 = std::sqrt(0.5);
    const double lvl1 = 1.0;
    const double thr  = 0.5 * (lvl0 + lvl1);

    QSBResult out;
    out.k2 = Eigen::MatrixXi::Zero(absD.rows(), absD.cols());
    out.p2 = Eigen::MatrixXd::Zero(absD.rows(), absD.cols());

    for (int i = 0; i < absD.rows(); ++i)
    {
        for (int t = 0; t < absD.cols(); ++t)
        {
            if (absD(i,t) >= thr)
            {
                out.k2(i,t) = 1;
                out.p2(i,t) = lvl1;
            }
            else
            {
                out.k2(i,t) = 0;
                out.p2(i,t) = lvl0;
            }
        }
    }

    return out;
}

struct QPResult
{
    Eigen::MatrixXi c;
    Eigen::MatrixXcd phi;
};

QPResult Q_P(const Eigen::MatrixXd& ang,
              int N_PSK)
{
    const int rows = ang.rows();
    const int SB   = ang.cols();

    QPResult out;
    out.c   = Eigen::MatrixXi(rows, SB);
    out.phi = Eigen::MatrixXcd(rows, SB);

    for (int i = 0; i < rows; ++i)
    {
        for (int t = 0; t < SB; ++t)
        {
            double a = std::fmod(ang(i,t), 2*M_PI);
            if (a < 0) a += 2*M_PI;

            int c = static_cast<int>(
                std::floor(a * N_PSK / (2*M_PI) + 0.5)
            ) % N_PSK;

            out.c(i,t) = c;

            double phase = 2*M_PI * c / N_PSK;
            out.phi(i,t) =
                std::exp(std::complex<double>(0.0, phase));
        }
    }

    return out;
}

Eigen::MatrixXcd BuildWC2(const Eigen::MatrixXcd& D,
                          int N_PSK,
                          Eigen::MatrixXd& P_SB)
{
    Eigen::MatrixXd absD = D.cwiseAbs();
    Eigen::MatrixXd angD = D.array().arg();

    auto qsb = Q_SB(absD);
    auto qp  = Q_P(angD, N_PSK);

    P_SB = qsb.p2;

    Eigen::MatrixXcd WC2 =
        qsb.p2.array() * qp.phi.array();

    return WC2;
}

static inline Eigen::DiagonalMatrix<std::complex<double>, Eigen::Dynamic>
BuildSigmaSB(const Eigen::VectorXd& p,          // (2L)
             const Eigen::MatrixXd& p2,         // (2L x SB)
             int N1, int N2)
{
    const int twoL = static_cast<int>(p.size());
    const int SB   = static_cast<int>(p2.cols());

    Eigen::VectorXcd diag(SB);
    const double scale = static_cast<double>(N1) * static_cast<double>(N2);

    for (int t = 0; t < SB; ++t)
    {
        double temp = 0.0;
        for (int j = 0; j < twoL; ++j)
        {
            // temp += p[j] * (p2[j,t]^2)
            const double a = p2(j, t);
            temp += p(j)*p(j) * (a * a);
        }

        // guard against zero/negative (shouldn't be negative)
        const double denom = std::max(scale * temp, 1e-30);
        diag(t) = 1.0 / std::sqrt(denom);
    }

    return Eigen::DiagonalMatrix<std::complex<double>, Eigen::Dynamic>(diag);
}

static inline Eigen::DiagonalMatrix<std::complex<double>, Eigen::Dynamic>
BuildWC1(const Eigen::VectorXd& p)
{
    Eigen::VectorXcd d(p.size());
    for (int i = 0; i < p.size(); ++i)
    {
        d(i) = std::complex<double>(p(i), 0.0);
    }
    return Eigen::DiagonalMatrix<std::complex<double>, Eigen::Dynamic>(d);
}

ComplexMatrixArray
BuildPrecMat(const Eigen::MatrixXcd& W,
             const Eigen::MatrixXcd& W_2,
             uint8_t m_rankMaxTp2)
{
    const int rows = W.rows();   // 2*Nant
    const int SB   = W.cols();

    NS_ASSERT_MSG(W_2.rows() == rows && W_2.cols() == SB, "W_2 size mismatch");

    const int rank = (m_rankMaxTp2 >= 2) ? 2 : 1;

    // PrecMat: rows x rank x SB
    ComplexMatrixArray precMat(rows, rank, SB);

    for (int t = 0; t < SB; ++t)
    {
        // Map the page t as (rows x rank)
        Eigen::Map<Eigen::MatrixXcd> page(precMat.GetPagePtr(t), rows, rank);

        // Column 0 always present
        page.col(0) = W.col(t);

        // Column 1 only if rank=2
        if (rank == 2)
        {
            page.col(1) = W_2.col(t);
        }
    }

    return precMat;
}

TypeId
NrPmSearchTypeTwo::GetTypeId()
{
    static TypeId tid = 
        TypeId("ns3::NrPmSearchTypeTwo")
            .SetParent<NrPmSearch>()
            .AddConstructor<NrPmSearchTypeTwo>()
            .AddAttribute("OversampH",
                        "Oversampling in the Horizontal port {Default to false}",
                        BooleanValue(true),
                        MakeBooleanAccessor(&NrPmSearchTypeTwo::m_oversampH),
                        MakeBooleanChecker())
            .AddAttribute("OversampV",
                        "Oversampling in the Vertical port {Default to false}",
                        BooleanValue(true),
                        MakeBooleanAccessor(&NrPmSearchTypeTwo::m_oversampV),
                        MakeBooleanChecker())
            .AddAttribute("BeamNumber",
                        "Number of beams in a layer for Type-II Codebook {1, 2, 3, 4}",
                        UintegerValue(2),
                        MakeUintegerAccessor(&NrPmSearchTypeTwo::m_numBeams),
                        MakeUintegerChecker<uint8_t>(1, 4));
    return tid;
}

void
NrPmSearchTypeTwo::InitCodebooks()
{
    auto maxRank = std::min({m_nRxPorts, m_nGnbPorts, static_cast<size_t>(m_rankLimit)});
    m_ranks.resize(maxRank);
    std::iota(m_ranks.begin(), m_ranks.end(), 1);
    m_rankParams.resize(maxRank + 1);
}

// ----------------------------------------------- MAIN Type-II Codebooks Implementation --------------------------------------------------------
PmCqiInfo
NrPmSearchTypeTwo::CreateCqiFeedbackMimo(const NrMimoSignal& rxSignalRb, PmiUpdate pmiUpdate)
{
    NS_LOG_FUNCTION(this);

    // Extract parameters from the received signal
    auto nRows = rxSignalRb.m_chanMat.GetNumRows();
    auto nCols = rxSignalRb.m_chanMat.GetNumCols();
    NS_ASSERT_MSG(nRows == m_nRxPorts, "Channel mat has {} rows but UE has {} ports");
    NS_ASSERT_MSG(nCols == m_nGnbPorts, "Channel mat has {} cols but gNB has {} ports");

    // Compute the interference-normalized channel matrix
    auto rbNormChanMat = rxSignalRb.m_covMat.CalcIntfNormChannel(rxSignalRb.m_chanMat);
    // Get the simulation time to dump the channel matrix into the CSV file
    const double t = Simulator::Now().GetSeconds();

    // DumpComplexMatrixArrayToCsv(
    // rbNormChanMat,
    // t,
    // "chanmat_dump.csv");


    // Dump the channel matrix into the file 
    /////// ---------- PLACEHOLDER ---------------

    // Update the optimal precoding matrices based on the received signal, if update is requested
    ConditionallyUpdatePrecoding(rbNormChanMat, pmiUpdate);

    auto OptPrecForRanks = std::vector<PmCqiInfo>{};
    for (auto rank : m_ranks)
    {
        auto cqiMsg = CreateCqiForRank(rank, rbNormChanMat);
        OptPrecForRanks.emplace_back(std::move(cqiMsg));

        // Skip higher ranks when the current is incapble of maintaining the connection
        if (OptPrecForRanks.back().m_wbCqi == 0)
        {
            if (OptPrecForRanks.size() >= 2)
            {
                OptPrecForRanks.pop_back();
            }
            break;
        }
    }

    // Find the rank which results in largest expected TB size, and return corresponding CQI/PMI
    auto optRankCqiMsg = *std::max_element(
        OptPrecForRanks.begin(),
        OptPrecForRanks.end(),
        [](const PmCqiInfo& a, const PmCqiInfo& b) {return a.m_tbSize < b.m_tbSize;}
    );

    return optRankCqiMsg;
}

PmCqiInfo
NrPmSearchTypeTwo::CreateCqiForRank(uint8_t rank, const NrIntfNormChanMat& rbNormChanMat) const
{
    auto optPrec = m_rankParams[rank].precParams;
    NS_ASSERT_MSG(optPrec, "Tried to create a CQI message but precoding matrix does not exist");

    // Upsample subband precoding matrix to full RB size
    auto rbPrecMat = SubbandUpsampling(optPrec->sbPrecMat, rbNormChanMat.GetNumPages());

    // Recompute the SINR value for the current channel (for all RBs)
    auto sinrMat = rbNormChanMat.ComputeSinrForPrecoding(rbPrecMat);

    // For optimal precoding matrix, determine the TB szie and TBLER.
    auto mcsParams = m_amc->GetMaxMcsParams(sinrMat, m_subbandSize);

    // Clamp sibband CQI according to 3GPP 2-bit overhead limit
    if (m_subbandCqiClamping)
    {
        for (auto& sbCqi : mcsParams.sbCqis)
        {
            auto diff = (int)sbCqi - (int)mcsParams.wbCqi;
            if (diff > 2)
            {
                sbCqi = mcsParams.wbCqi + 2;
            }
            else if (diff < -1)
            {
                sbCqi = mcsParams.wbCqi - 1;
            }
        }
    }

    // Store parameters for this rank
    // auto cqiMsg = PmCqiInfo{
    //     .m_mcs = mcsParams.mcs,
    //     .m_rank = rank,
    //     .m_wbCqi = mcsParams.wbCqi,
    //     .m_sbCqis = mcsParams.sbCqis,
    //     .m_optPrecMat = Create<const ComplexMatrixArray>(std::move(rbPrecMat)),
    //     .m_tbSize = mcsParams.tbSize,
    // };

    auto cqiMsg = PmCqiInfo{
        .m_mcs = mcsParams.mcs,
        .m_rank = rank,
        .m_wbPmi = optPrec->wbPmi,
        .m_wbCqi = mcsParams.wbCqi,
        .m_sbCqis = mcsParams.sbCqis,
        .m_sbPmis = optPrec->sbPmis,
        .m_optPrecMat = Create<const ComplexMatrixArray>(std::move(rbPrecMat)),
        .m_tbSize = mcsParams.tbSize,
    };

    return cqiMsg;
}

void
NrPmSearchTypeTwo::ConditionallyUpdatePrecoding(const NrIntfNormChanMat& rbNormChanMat,
                                                PmiUpdate pmiUpdate)
{
    if (pmiUpdate.updateWb || pmiUpdate.updateSb)
    {
        UpdateAllPrecoding(rbNormChanMat);
    }
    // else if (pmiUpdate.updateSb)
    // {
    //     // UpdateSubbandPrecoding(rbNormChanMat);
    //     UpdateAllPrecoding(rbNormChanMat);
    // }
}

void
NrPmSearchTypeTwo::UpdateAllPrecoding(const NrIntfNormChanMat& rbNormChanMat)
{
    // Compute downsampled channel per subband
    auto sbNormChanMat = SubbandDownsampling(rbNormChanMat);

    // Step 1: Compute the long-term precoding matrix
    // Compute the matrix R of size (nGnbPorts x n GnbPorts x Subbands)
    auto Hcorr = NrIntfNormChanMat(sbNormChanMat.HermitianTranspose() * sbNormChanMat);

    // Select the maximum rank
    // m_maxRankTp2 = SelectRank(Hcorr);

    // Compue the average correlation matrix R of the channel matrix of size (nGnbPorts x nGnbPorts)
    size_t n = Hcorr.GetNumRows();
    size_t nPages = Hcorr.GetNumPages();
    // auto R = ComplexMatrixArray{n, n};

    Eigen::MatrixXcd Havg = Eigen::MatrixXcd::Zero(n, n);

    for (size_t rb = 0; rb < nPages; ++rb)
    {
        Eigen::Map<const Eigen::MatrixXcd> Hrb(
            reinterpret_cast<const std::complex<double>*>(Hcorr.GetPagePtr(rb)),
            n,
            n
        );

        Havg += Hrb;
    }

    Havg /= static_cast<double>(nPages);
    // std::cout << "Channel matrix correlation matrix dimension: HAVG" << std::endl;
    // std::cout<< "ROWS = " << Havg.rows() << std::endl;
    // std::cout<< "COLS = " << Havg.cols() << std::endl;

    // Divide the block matrix and extract the main diag parts
    size_t N = n / 2;
    Eigen::MatrixXcd R11 = Havg.block(0, 0, N, N);
    Eigen::MatrixXcd R12 = Havg.block(0, N, N, N);
    Eigen::MatrixXcd R21 = Havg.block(N, 0, N, N);
    Eigen::MatrixXcd R22 = Havg.block(N, N, N, N);

    NS_ASSERT_MSG((R21 - R12.adjoint()).norm() < 1e-9, "The cross diag should be orthogonal");

    Eigen::MatrixXcd R = R11 + R22;


    // Find the best long-term beams based on the channel correlation matrix
    uint8_t O1 = m_oversampH ? 4 : 1;
    uint8_t O2 = m_oversampV ? 4 : 1;

    std::cout << "The H and V ports = " << m_nGnbHPorts << ", " << m_nGnbVPorts << std::endl;
    
    auto res = SelectBeamsTypeTwo(R,
                                    /*N1*/ m_nGnbHPorts, /*N2*/ m_nGnbVPorts,
                                    /*O1*/ O1, /*O2*/ O2,
                                    /*L*/ m_numBeams );


    // Build the W1 matrix
    const uint32_t L = res.v.size();
    const uint32_t Nant = res.v[0].rows();

    Eigen::MatrixXcd W1 = Eigen::MatrixXcd::Zero(2 * Nant, 2 * L);
    for (uint32_t i = 0; i < L; ++i)
    {
        W1.block(0,        i, Nant, 1) = res.v[i];
        W1.block(Nant, L + i, Nant, 1) = res.v[i];
    }

    // Step 2: Build the effective channel
    ComplexMatrixArray HA = RightMultiply(sbNormChanMat, W1);

    // Step 3: Compute the Normalize average W_2 matrix
    auto resW = ComputeW2andNormalize(HA);

    // Step 4: Compute wideband and subband amplitudes and phase
    for (auto rank : m_ranks)
    {
        if (rank == 1)
        {
            auto precRes = Create<PrecMatParams>();
            Eigen::VectorXd m(resW.W2.rows());

            for (int r = 0; r < resW.W2.rows(); ++r)
            {
                m(r) = resW.W2.row(r).cwiseAbs().maxCoeff();
            }

            auto qres = Q_WB(m);

            // Eigen::VectorXcd k = qres.k;
            Eigen::VectorXd p = qres.p;

            Eigen::MatrixXcd D = ComputeD(resW.W2, p);

            Eigen::MatrixXd P_SB;
            Eigen::MatrixXcd WC2 = BuildWC2(D, 4, P_SB);

            auto Sigma_SB = BuildSigmaSB(p, P_SB, m_nGnbHPorts, m_nGnbVPorts);
            auto WC1  = BuildWC1(p);

            Eigen::MatrixXcd W  = W1 * WC1  * WC2  * Sigma_SB;

            ComplexMatrixArray PrecMat = BuildPrecMat(W, W, rank);

            precRes->sbPrecMat = PrecMat;

            m_rankParams[rank].precParams = precRes;
        }
        else if (rank == 2)
        {
            auto precRes = Create<PrecMatParams>();
            Eigen::VectorXd m(resW.W2.rows());
            Eigen::VectorXd m2(resW.W22.rows());

            for (int r = 0; r < resW.W2.rows(); ++r)
            {
                m(r) = resW.W2.row(r).cwiseAbs().maxCoeff();
                m2(r) = resW.W22.row(r).cwiseAbs().maxCoeff();
            }

            auto qres = Q_WB(m);
            auto qres2 = Q_WB(m2);

            // Eigen::VectorXcd k = qres.k;
            // Eigen::VectorXcd k2 = qres2.k;

            Eigen::VectorXd p = qres.p;
            Eigen::VectorXd p2 = qres2.p;

            Eigen::MatrixXcd D = ComputeD(resW.W2, p);
            Eigen::MatrixXcd D2 = ComputeD(resW.W22, p2);

            Eigen::MatrixXd P_SB;
            Eigen::MatrixXd P_SB2;

            Eigen::MatrixXcd WC2 = BuildWC2(D, 4, P_SB);
            Eigen::MatrixXcd WC22 = BuildWC2(D2, 4, P_SB2);

            auto Sigma_SB = BuildSigmaSB(p, P_SB, m_nGnbHPorts, m_nGnbVPorts);
            auto WC1  = BuildWC1(p);

            auto Sigma_SB2 = BuildSigmaSB(p2, P_SB2, m_nGnbHPorts, m_nGnbVPorts);
            auto WC12  = BuildWC1(p2);

            Eigen::MatrixXcd W  = W1 * WC1  * WC2  * Sigma_SB;
            Eigen::MatrixXcd W_2  = W1 * WC12  * WC22  * Sigma_SB2;

            ComplexMatrixArray PrecMat2 = BuildPrecMat(W, W_2, rank);

            precRes->sbPrecMat = PrecMat2;

            m_rankParams[rank].precParams = precRes;
        }
        
    }

}

} // namespace ns3