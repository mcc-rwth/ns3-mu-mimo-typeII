#ifndef NR_PM_SEARCH_TYPE_TWO_H
#define NR_PM_SEARCH_TYPE_TWO_H

#include "nr-amc.h"
#include "nr-mimo-chunk-processor.h"
#include "nr-mimo-signal.h"

#include "ns3/random-variable-stream.h"
#include "ns3/uinteger.h"

#include "nr-cb-type-two.h"
#include "ns3/object-factory.h"

#include "nr-pm-search.h"

namespace ns3
{

/**
 * @brief Class for searching optimal precoding matrices and creating full CQI/PMI feedback
 */

 class NrPmSearchTypeTwo : public NrPmSearch
 {
    public:
        /**
         * @brief Get TypeId
         * @return The TypeId
         */
        static TypeId GetTypeId();

        /**
         *  @brief Create and initialize the codebook for each rank.
         */
        virtual void InitCodebooks() override;

        /**
         * @brief Create CQI feedback with optimal rank, optimal PMI, and corresponding CQI values.
         * Optimal rank is considered as the rank tham maximizes the achievable TB size when using the
         * optimal PMI. The optimal WB/SB PMI values are updated based on PmiUpdate. If there is no update
         * to the PMI values, the previously found PMI values can be used.
         * @param rxSignalRb The receive signal parameters (channel and interference matrices)
         * @param pmiUpdate struct that defines if WB/SB PMIs need to be updated
         * @return The CQI feedback message that contains the optimum RI, PMI, CQI, and full precoding matrix
         * Dimensions: nGnbPorts x rank x nRbs
         */
        PmCqiInfo CreateCqiFeedbackMimo(const NrMimoSignal& rxSignalRb,
                                                PmiUpdate PmiUpdate) override;

    protected: 
        struct RankParams
        {
            Ptr<PrecMatParams> precParams;       ///< The precoding parameters (WB/SB PMIs)
            // Ptr<NrCbTypeTwo> cb;                    ///< The Codebook
        };



        /**
         * @brief Update the WB and/or SB PMI, or neither.
         * @param rbNormChanMat the interference-normed channel matrix per RB
         * @param pmiUpdate the struct defining if updates to SB or WB PMI are necessary
         */
        void ConditionallyUpdatePrecoding(const NrIntfNormChanMat& rbNormChanMat, PmiUpdate pmiUpdate);

        /**
         * @brief For all ranks, update the optimum precoding matrices (wideband and subband).
         * @param rbNormChanMat the interference-normed channel matrix per RB
         */
        void UpdateAllPrecoding(const NrIntfNormChanMat& rbNormChanMat);

        // /**
        //  * @brief For all ranks, update the opt subband PMI assuming previous value of wideband PMI.
        //  * @param rbNormChanMat the interference-normed channel matrix per RB
        //  */
        // void UpdateSubbandPrecoding(const NrIntfNormChanMat& rbNormChanMat);

        /**
         * @brief Create CQI feedback message for a particular rank.
         * @param rank the rank for which to create the CQI feedback
         * @param rbNormChanMat the interference-normed channel matrix per RB
         * @return the CQI message with PMI and CQI values, as well as expected TB size
         */
        PmCqiInfo CreateCqiForRank(uint8_t rank, const NrIntfNormChanMat& rbNormChanMat) const;
        
        
        bool m_oversampV{false};      ///< Antenna vertical oversampling
        bool m_oversampH{false};      ///< Antenna horizontal oversampling
        size_t m_numBeams{1};       ///< Number of beams in a layer

        std::vector<RankParams> m_rankParams;   ///< The parameters (PMI values, codebook) for each rank
        ObjectFactory m_cbfactory;              ///< The factory used to create the codebooks
        uint8_t m_maxRankTp2;                   ///< Selected maximum rank of the precoding based on the channel matrix
    };
} // namespace ns3

#endif // NR_PM_SEARCH_TYPE_TWO_H