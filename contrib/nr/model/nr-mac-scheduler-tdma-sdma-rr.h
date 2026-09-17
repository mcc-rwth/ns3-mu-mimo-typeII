#pragma once

#include "nr-mac-scheduler-tdma-sdma.h"

#include <unordered_set>

namespace ns3
{

class NrMacSchedulerTdmaSdmaRR : public NrMacSchedulerTdmaSdma
{
    public:

        static TypeId GetTypeId();

        NrMacSchedulerTdmaSdmaRR();

        ~NrMacSchedulerTdmaSdmaRR() override
        {
        }

    protected:

        std::shared_ptr<NrMacSchedulerUeInfo> CreateUeRepresentation(
                const NrMacCschedSapProvider::CschedUeConfigReqParameters& params) const override;

        /**
         * @brief Return the comparison function to sort DL UE according to the scheduler policy
         * @return a pointer to NrMacSchedulerUeInfoRR::CompareUeWeightsDl
         */
        std::function<bool(const NrMacSchedulerNs3::UePtrAndBufferReq& lhs,
                        const NrMacSchedulerNs3::UePtrAndBufferReq& rhs)>
        GetUeCompareDlFn() const override;

        /**
         * @brief Return the comparison function to sort UL UE according to the scheduler policy
         * @return a pointer to NrMacSchedulerUeInfoRR::CompareUeWeightsUl
         */
        std::function<bool(const NrMacSchedulerNs3::UePtrAndBufferReq& lhs,
                        const NrMacSchedulerNs3::UePtrAndBufferReq& rhs)>
        GetUeCompareUlFn() const override;

        /**
         * @brief Update the UE representation after a symbol (DL) has been assigned to it
         * @param ue UE to which a symbol has been assigned
         * @param assigned the amount of resources assigned
         * @param totAssigned the total amount of resources assigned in the slot
         *
         * Update DL metrics by calling NrMacSchedulerUeInfoRR::UpdateDlMetric
         */
        void AssignedDlResources(const UePtrAndBufferReq& ue,
                                const FTResources& assigned,
                                const FTResources& totAssigned) const override;

        /**
         * @brief Update the UE representation after a symbol (DL) has been assigned to it
         * @param ue UE to which a symbol has been assigned
         * @param assigned the amount of resources assigned
         * @param totAssigned the total amount of resources assigned in the slot
         *
         * Update DL metrics by calling NrMacSchedulerUeInfoRR::UpdateUlMetric
         */
        void AssignedUlResources(const UePtrAndBufferReq& ue,
                                const FTResources& assigned,
                                const FTResources& totAssigned) const override;

        // RR is a simple scheduler: it doesn't do anything in the next
        // inherited calls.
        void NotAssignedDlResources(const UePtrAndBufferReq& ue,
                                    const FTResources& notAssigned,
                                    const FTResources& totalAssigned) const override
        {
        }

        void NotAssignedUlResources(const UePtrAndBufferReq& ue,
                                    const FTResources& notAssigned,
                                    const FTResources& totalAssigned) const override
        {
        }

        void BeforeDlSched(const UePtrAndBufferReq& ue,
                        const FTResources& assignableInIteration) const override
        {
        }

        void BeforeUlSched(const UePtrAndBufferReq& ue,
                        const FTResources& assignableInIteration) const override
        {
        }

    private:
        mutable std::deque<uint16_t> m_dlRrRntiDeque;
        mutable std::unordered_set<uint16_t> m_dlRntiSet;

};
} // namespace ns3