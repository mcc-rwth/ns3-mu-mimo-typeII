#include "nr-cb-type-two.h"

#include "ns3/double.h"
#include "ns3/enum.h"
#include "ns3/log.h"
#include "ns3/math.h"
#include "ns3/nr-spectrum-value-helper.h"
#include "ns3/uinteger.h"

#include "complex.h"

namespace ns3
{

    using namespace std::complex_literals;

    NS_LOG_COMPONENT_DEFINE("NrCbTypeTwo");
    NS_OBJECT_ENSURE_REGISTERED(NrCbTypeTwo);

    TypeId
    NrCbTypeTwo::GetTypeId()
    {
        static TypeId tid = 
                    TypeId("ns3::NrCbTypeTwo").SetParent<NrCbTypeTwoBase>().AddConstructor<NrCbTypeTwo>();
        
        return tid;
    }

    void
    NrCbTypeTwo::Init()
    {
        NS_ASSERT_MSG(m_rank > 0, "Rank must not be zero!");
        NS_ASSERT_MSG(m_rank <= 2, "The codebook supports at most rank 2 MIMO!");

        m_o1 = (m_n1 > 0) ? 4 : 1;
        m_o2 = (m_n2 > 0) ? 4 : 1;
        m_nPorts = (m_isDualPol) ? 2 * m_n1 * m_n2 : m_n1 * m_n2;

        NS_ASSERT_MSG(m_nPorts > 0, "Number of CSI-RS ports must not be zero!");
        NS_ASSERT_MSG(m_isDualPol || m_nPorts <= 2, "For > 2 antenna ports, dual polarization is required.");
        NS_ASSERT_MSG(m_nPorts <= 32, "Number of CSI-RS ports must not be greater than 32");

        InitWParams();
    }

    ComplexMatrixArray
    NrCbTypeTwo::GetBasePrecMat(std::vector<size_t> i1, std::vector<std::vector<size_t>> i2) const
    {
        auto i11 = MapToI11(i1);
        auto i12 = MapToI12(i1);
        auto i13 = MapToI13(i1);
        auto i14 = MapToI14(i1);

        auto i21 = MapToI21(i2);
        auto i22 = MapToI22(i2);

        return GetBasePrecMatFromIndex(i11, i12, i13, i14, i21, i22);
    }

    ComplexMatrixArray
    NrCbTypeTwo::GetBasePrecMatFromIndex(size_t i11, size_t i12, size_t i13, std::vector<size_t> i14, size_t i21, size_t i22) const
    {
        if (m_nPorts == 1)
        {
            auto res = ComplexMatrixArray(1, 1);
            res(0, 0) = 1.0;
            return res;
        }

        // m_nPorts is even-numbered. The upper half of the ports represent the first polarization angle
        auto precMat = ComplexMatrixArray(m_nPorts, m_rank);
        // auto phase = M_PI * static_cast<double> Fix this for Type-II codebook case of the co-phasing values!
        auto normalizer = 1.0 / sqrt(m_nPorts * m_rank); // TODO: Fix this for Type-II codebook case
        auto uniqueBfvs = CreateUniqueBfvs(i11, i12, i13, i14);

        for (size_t stream = 0; stream < m_rank; stream++)
        {
            // The beamforming vector for the first polarization
            // auto v = uniqueBfvs[];
        }
    }
    
    size_t
    NrCbTypeTwo::Cfn(size_t x, size_t y) const
    {
        // Lower triangular always zero
        if (y > x) return 0;

        if (y == 0 || x == y) return 1;

        double res = 0;
        for (int i = 0; i < y; i++)
        {
            res += std::log(x - i) - std::log(i + 1);
        }

        return (int)round(std::exp(res));
    }
}