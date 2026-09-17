#include "nr-cb-type-two-base.h"

#include "ns3/boolean.h"
#include "ns3/uinteger.h"

namespace ns3
{

NS_LOG_COMPONENT_DEFINE("NrCbTypeTwoBase");
NS_OBJECT_ENSURE_REGISTERED(NrCbTypeTwoBase);

TypeId
NrCbTypeTwoBase::GetTypeId()
{
    static TypeId tid = TypeId("ns3::NrCbTypeTwoBase")
                            .SetParent<Object>()
                            .AddAttribute("N1",
                                        "N1 (the number of horizontal ports at the gNB)",
                                        UintegerValue(NR_CB_TYPE_TWO_INIT_N1),
                                        MakeUintegerAccessor(&NrCbTypeTwoBase::m_n1),
                                        MakeUintegerChecker<uint8_t>())
                            .AddAttribute("N2",
                                        "N2 (the number of vertical ports at the gNB)",
                                        UintegerValue(NR_CB_TYPE_TWO_INIT_N2),
                                        MakeUintegerAccessor(&NrCbTypeTwoBase::m_n2),
                                        MakeUintegerChecker<uint8_t>())
                            .AddAttribute("IsDualPol",
                                        "True if the gNB antennas are dual polarized",
                                        BooleanValue(NR_CB_TYPE_TWO_INIT_POL),
                                        MakeBooleanAccessor(&NrCbTypeTwoBase::m_isDualPol),
                                        MakeBooleanChecker())
                            .AddAttribute("Rank",
                                        "Rank (the number of MIMO transmissions)",
                                        UintegerValue(NR_CB_TYPE_TWO_INIT_RANK),
                                        MakeUintegerAccessor(&NrCbTypeTwoBase::m_rank),
                                        MakeUintegerChecker<uint8_t>(1, 2));
    return tid;
}

size_t
NrCbTypeTwoBase::GetNumI1() const
{
    return m_numI1;
}

size_t
NrCbTypeTwoBase::GetNumI2() const
{
    return m_numI2;
}

} // namespace ns3