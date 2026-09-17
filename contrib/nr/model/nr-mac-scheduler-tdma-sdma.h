#pragma once

#include "nr-mac-scheduler-ns3.h"

#include <functional>
#include <memory>
#include <unordered_set>

namespace ns3
{
/**
 * @ingroup scheduler
 * @brief The SDMA-TDMA scheduler for MU-MIMO
 * 
 * An example of SDMA-TDMA-based scheduling is the following
 * <pre>
 * (f)
 * ^
 * |=|===|===|===|===|=|
 * |C| U | U | U | U |C|
 * |T| E | E | U | U |T|
 * |R| G | G | G | G |R|
 * |L| 1 | 2 | 3 | 4 |L|
 * |--------------------------------> (t)
 * <pre>
 */

class NrMacSchedulerTdmaSdma : public NrMacSchedulerNs3
{
    friend class NrTestSchedulerAiCase;

    public:
        /**
         * @brief GetTypeId
         * @return the TypeId of the class
         */
        static TypeId GetTypeId();

        NrMacSchedulerTdmaSdma();

        ~NrMacSchedulerTdmaSdma() override;

        enum GroupingMetric : uint8_t
        {
            GROUP_ZF_DIR = 0,
            GROUP_COLINEAR = 1,
            GROUP_VOLUME = 2,
            GROUP_KMEANS = 3,
        };


    protected:

        BeamSymbolMap AssignDLRBG(uint32_t symAvail, const ActiveUeMap& activeDl) const override;
        BeamSymbolMap AssignULRBG(uint32_t symAvail, const ActiveUeMap& activeUl) const override;

        std::shared_ptr<DciInfoElementTdma> CreateDlDci(
            PointInFTPlane* spoint,
            const std::shared_ptr<NrMacSchedulerUeInfo>& ueInfo,
            uint32_t maxSym) const override;

        std::shared_ptr<DciInfoElementTdma> CreateDlDci(
            PointInFTPlane* spoint,
            const std::shared_ptr<NrMacSchedulerUeInfo>& ueInfo,
            uint32_t maxSym, bool advanceSpoint) const;

        std::shared_ptr<DciInfoElementTdma> CreateUlDci(
            PointInFTPlane* spoint,
            const std::shared_ptr<NrMacSchedulerUeInfo>& ueInfo,
            uint32_t maxSym) const override;

        /**
         * @brief Not doing anything, moving forward the spoint is done by CreateDci
         * @param spoint Starting point
         * @param symOfBeam the number of symbols assigned to the beam
         */
        void ChangeDlBeam([[maybe_unused]] PointInFTPlane* spoint,
                        [[maybe_unused]] uint32_t symOfBeam) const override
        {
        }

        /**
         * @brief Not doing anything, moving forward the spoint is done by CreateDci
         * @param spoint Starting point
         * @param symOfBeam the number of symbols assigned to the beam
         */
        void ChangeUlBeam([[maybe_unused]] PointInFTPlane* spoint,
                        [[maybe_unused]] uint32_t symOfBeam) const override
        {
        }

        /**
         * @brief Provide the comparison function to order the UE when scheduling DL
         * @return a function that should order two UEs based on their priority: if
         * UE a is less than UE b, it will have an higher priority.
         */
        virtual std::function<bool(const NrMacSchedulerNs3::UePtrAndBufferReq& lhs,
                                const NrMacSchedulerNs3::UePtrAndBufferReq& rhs)>
        GetUeCompareDlFn() const = 0;

        /**
         * @brief Provide the comparison function to order the UE when scheduling UL
         * @return a function that should order two UEs based on their priority: if
         * UE a is less than UE b, it will have an higher priority.
         */
        virtual std::function<bool(const NrMacSchedulerNs3::UePtrAndBufferReq& lhs,
                                const NrMacSchedulerNs3::UePtrAndBufferReq& rhs)>
        GetUeCompareUlFn() const = 0;

        virtual void AssignedDlResources(const UePtrAndBufferReq& ue,
                                            const FTResources& assigned,
                                            const FTResources& totalAssigned) const = 0;

        /**
         * @brief Update the UE representation after a symbol (DL) has been assigned to it
         * @param ue UE to which a symbol has been assigned
         * @param assigned the amount of resources assigned
         * @param totalAssigned the amount of total resources assigned until now
         *
         * After an UE is selected to be eligible for a symbol assignment, its representation
         * should be updated. The subclasses, by implementing this method, update
         * the representation by updating some custom values that reflect the assignment
         * done. These values are the one that, hopefully, are checked by the
         * comparison function returned by GetUeCompareUelFn().
         */
        virtual void AssignedUlResources(const UePtrAndBufferReq& ue,
                                        const FTResources& assigned,
                                        const FTResources& totalAssigned) const = 0;

        /**
         * @brief Update the UE representation after a symbol (DL) has been assigned to other UE
         * @param ue UE to which a symbol has not been assigned
         * @param notAssigned the amount of resources not assigned
         * @param totalAssigned the amount of total resources assigned until now
         */
        virtual void NotAssignedDlResources(const UePtrAndBufferReq& ue,
                                            const FTResources& notAssigned,
                                            const FTResources& totalAssigned) const = 0;


        /**
         * @brief Update the UE representation after a symbol (UL) has been assigned to other UE
         * @param ue UE to which a symbol has not been assigned
         * @param notAssigned the amount of resources not assigned
         * @param totalAssigned the amount of total resources assigned until now
         */
        virtual void NotAssignedUlResources(const UePtrAndBufferReq& ue,
                                            const FTResources& notAssigned,
                                            const FTResources& totalAssigned) const = 0;

        /**
        * @brief Prepare UE for the DL scheduling
        * @param ue UE that is eligible for an assignation in any iteration round
        * @param assignableInIteration Resources that can be assigned in each iteration
        *
        * The default implementation is empty, but a subclass can specialize the
        * behaviour, e.g., to calculate some value before the choice of RBG to
        * assign to each UE is done.
        */
        virtual void BeforeDlSched(const UePtrAndBufferReq& ue,
                                const FTResources& assignableInIteration) const = 0;

        /**
        * @brief Prepare UE for the UL scheduling
        * @param ue UE that is eligible for an assignation in any iteration round
        * @param assignableInIteration Resources that can be assigned in each iteration
        *
        * The default implementation is empty, but a subclass can specialize the
        * behaviour, e.g., to calculate some value before the choice of RBG to
        * assign to each UE is done.
        */
        virtual void BeforeUlSched(const UePtrAndBufferReq& ue,
                                const FTResources& assignableInIteration) const = 0;

        /**
         * @brief Call the notify callback function in the OpenGymEnv class
         * in the ns3-gym module for downlink
         * @param ueVector A vector containing pointers to active UEs and their corresponding buffer
         * requests
         */
        virtual void CallNotifyDlFn(const std::vector<UePtrAndBufferReq>& ueVector) const
        {
        }

        /**
         * @brief Call the notify callback function in the OpenGymEnv class
         * in the ns3-gym module for uplink
         * @param ueVector A vector containing pointers to active UEs and their corresponding buffer
         * requests
         */
        virtual void CallNotifyUlFn(const std::vector<UePtrAndBufferReq>& ueVector) const
        {
        }

        typedef std::function<bool(const NrMacSchedulerNs3::UePtrAndBufferReq& lhs,
                                    const NrMacSchedulerNs3::UePtrAndBufferReq& rhs)>
            CompareUeFn;
        typedef std::function<CompareUeFn()> GetCompareUeFn;


        /**
         * @brief Sorts the given vector of UE pointers and buffer requests.
         *
         * This method sorts the provided `ueVector` in-place using a stable sort algorithm, based on
         * the comparison function specified by `GetCompareFn`.
         *
         * Note: This function is inherited and overridden by nr-mac-scheduler-tdma-random to implement
         * random scheduler in which GetCompareFn is not used.
         *
         * @param ueVector A pointer to the vector containing UE pointers and their corresponding buffer
         * requests. The elements within this vector will be sorted according to the comparison criteria
         * defined by the `GetCompareFn`.
         * @param GetCompareFn A constant reference to a comparison function object that defines the
         * sorting order for the elements in the `ueVector`. This function should take two arguments of
         * type `UePtrAndBufferReq` and return a boolean indicating whether the first argument is less
         * than the second.
         */
        virtual void SortUeVector(std::vector<UePtrAndBufferReq>* ueVector,
                                [[maybe_unused]] const GetCompareUeFn& GetCompareFn) const;

        uint8_t GetTpc() const override;

        /**
         * @brief the function to check wether two UEs are compatible or not based on their precoding matrix W_i
         * @return bool return of compatibility check
         */
        bool IsCompatibleDl(const UePtrAndBufferReq& ueA, const UePtrAndBufferReq& ueB) const;

        double CompatibilityScore() const;

        void ComputeDlGroupPrecoding(const std::vector<UePtr>& group) const override;

        void ComputeFullCsiSvdPrecoding (const std::vector<UePtr>& g, uint32_t maxStreams) const;
        void ComputeAggregateSvdPrecoding(const std::vector<UePtr>& g, uint32_t maxStreams) const;
        // void FinalizeAllocation(std::vector<struct ns3::UeAlloc>& allocs, uint32_t Ntx, size_t pages) const;
        // ^ if your compiler dislikes the forward-declared struct, move UeAlloc
        //   out of the anonymous namespace into the .h, or make these free
          //   functions in the .cc instead of members.
        Ptr<const ComplexMatrixArray> GetDlChannelMatrix(const std::shared_ptr<NrMacSchedulerUeInfo>&) const;
        // uint8_t SeToMcs(double sePerLayer) const;

        
    private:
        /**
         * @brief Retrieve the UE vector from an ActiveUeMap
         * @param activeUes UE map
         * @return A Vector of UEs and their buffer requirements (in B)
         *
         * Really used only in TDMA scheduling. Worth moving?
         */
        static std::vector<UePtrAndBufferReq> GetUeVectorFromActiveUeMap(const ActiveUeMap& activeUes);


        typedef std::function<void(const UePtrAndBufferReq&, const FTResources&)>
            BeforeSchedFn; //!< Before scheduling function
        /**
         * @brief //!< Function to notify a successful assignment
         */
        typedef std::function<void(const UePtrAndBufferReq&, const FTResources&, const FTResources&)>
            AfterSuccessfulAssignmentFn;
        /**
         * @brief Function to notify that the UE did not get any resource in one iteration
         */
        typedef std::function<void(const UePtrAndBufferReq&, const FTResources&, const FTResources&)>
            AfterUnsuccessfulAssignmentFn;
        typedef std::function<std::vector<uint16_t>&(const UePtr& ue)>
            GetRBGFn;                                               //!< Getter for the RBG of an UE
        typedef std::function<uint32_t&(const UePtr& ue)> GetTBSFn; //!< Getter for the TBS of an UE
        typedef std::function<std::vector<uint8_t>&(const UePtr& ue)>
            GetSymFn; //!< Getter for the number of symbols of an UE
        typedef std::function<void(std::vector<UePtrAndBufferReq>& ueVector)> CallNotifyFn;

        BeamSymbolMap AssignRBGTDMA(uint32_t symAvail,
                                    const ActiveUeMap& activeUe,
                                    const std::string& type,
                                    const BeforeSchedFn& BeforeSchedFn,
                                    const GetCompareUeFn& GetCompareFn,
                                    const GetTBSFn& GetTBSFn,
                                    const GetRBGFn& GetRBGFn,
                                    const GetSymFn& GetSymFn,
                                    const AfterSuccessfulAssignmentFn& SuccessfulAssignmentFn,
                                    const AfterUnsuccessfulAssignmentFn& UnSuccessfulAssignmentFn,
                                    const CallNotifyFn& callNotifyFn) const;

        std::shared_ptr<DciInfoElementTdma> CreateDci(
            PointInFTPlane* spoint,
            const std::shared_ptr<NrMacSchedulerUeInfo>& ueInfo,
            uint32_t tbs,
            DciInfoElementTdma::DciFormat fmt,
            uint32_t mcs,
            uint8_t rank,
            Ptr<const ComplexMatrixArray> precMats,
            uint8_t numSym) const;

        std::vector<DciInfoElementTdma> DoReshapeAllocation(
            const std::vector<DciInfoElementTdma>& dcis,
            uint8_t& startingSymbol,
            uint8_t& numSymbols,
            std::vector<bool>& bitmask,
            const bool isDl,
            const std::unordered_map<uint16_t, std::shared_ptr<NrMacSchedulerUeInfo>>& ueMap) override;

        
        std::vector<std::vector<std::shared_ptr<NrMacSchedulerUeInfo>>>
            GroupByCollinearity(const std::vector<std::shared_ptr<NrMacSchedulerUeInfo>>& ues, uint32_t maxStreams) const;

        std::vector<std::vector<std::shared_ptr<NrMacSchedulerUeInfo>>>
            GroupByVolume(const std::vector<std::shared_ptr<NrMacSchedulerUeInfo>>& ues, uint32_t maxStreams) const;

        std::vector<std::vector<std::shared_ptr<NrMacSchedulerUeInfo>>>
            GroupByKMeans(const std::vector<std::shared_ptr<NrMacSchedulerUeInfo>>& ues, uint32_t maxStreams) const;

        std::vector<std::vector<std::shared_ptr<NrMacSchedulerUeInfo>>>
            GroupUesForSdma(const std::vector<std::shared_ptr<NrMacSchedulerUeInfo>>& ues) const;


        uint32_t m_maxSdmaUsers {4};
        double m_compatThresh {0.1}; ///!< e.g. |W_a^H W_b|^2 < thresh
        
    // =============================================================================
        bool     m_fullCsiBaseline {false};
        bool     m_aggregateSvd    {false};   // false = ZF, true = aggregate-SVD
        uint32_t m_maxStreams      {0};       // 0 = unconstrained; 4 = Part 2
        double   m_dlNoisePowerLin {1.0};     // (D) noise reference for MCS recompute

        uint32_t m_maxRankPerUe {1};
        double   m_sinrBackoffDb {3.0};
        int      m_forceDlMcs    {-1};
    
    // ======================================================================

        uint8_t m_groupingMetric {GROUP_ZF_DIR};
        double m_collinearityThr{0.30};
        double m_volumeThr{0.5};
        bool m_volumeUseAmGm{false};

    // ===============================================================================

        uint32_t m_kmeansClusters{0};
        uint32_t m_kmeansMaxIter{20};

};
} // namespace ns3