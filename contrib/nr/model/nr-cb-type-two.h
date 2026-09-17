#ifndef NR_CB_TYPE_TWO_H
#define NR_CB_TYPE_TWO_H

#include "nr-cb-type-two-base.h"

namespace ns3
{
class NrPmSearchFull;   ///< Forward declaration of friendship

/**
 * @brief Implementation of Type-II Codebook 3GPP TS 38.214 Rel. 18 Sec. 5.2.2.2.3
 * The rank is limited to 2.
 */
class NrCbTypeTwo : public NrCbTypeTwoBase
{
    public:
        /**
         * @brief Get ns-3 TypeId
         * @return the TypeId
         */
        static TypeId GetTypeId();

        /**
         * @brief Initialize the codebook parameters after construction, based on attributes values.
         */
        void Init() override;

        /**
         * @brief Get the 2D precoding matrix
         * @param i1 the composite index of the wideband precoding
         * @param i2 the composite index of the subband precoding
         * @return the precoding matrix of size m_nPorts x m_rank
         */
        ComplexMatrixArray GetBasePrecMat(std::vector<size_t> i1, std::vector<std::vector<size_t>> i2) const override;

        /**
         * @brief Get the 2D precoding matrix.
         * @param i11 the oversampling choice [q_1, q_2]
         * @param i12 the beam choice [n_1, n_2]
         * @param i13 the best beam for that layer
         * @param i14 wideband amplitude values
         * @param i21 subband phase values
         * @param i22 subband amplitude values
         */
        ComplexMatrixArray GetBasePrecMatFromIndex(size_t i11,
                                                size_t i12,
                                                size_t i13,
                                                std::vector<size_t> i14,
                                                size_t i21,
                                                size_t i22) const;

    
    protected:

        /// @brief Init the number of i11 indices (oversampling choice)
        void InitNumI11();

        /// @brief Init the number of i12 indices (beam choice)
        void InitNumI12();

        /// @brief Init the number of i13 indices (maximum beam)
        void InitNumI13();

        /// @brief Init the number of i14 indices (wideband amplitude)
        void InitNumI14();

        /**
         * @brief Define the columns in the full precoding matrix
         * @note Specifically, define for each column/layer the index of the vector to be selected from
         * uniqueBfvs, and the sign for the second polarization (sign in the lower half before phi_n)
         */
        void InitWParams();

        size_t MapToI11(std::vector<size_t> i1) const;
        
        size_t MapToI12(std::vector<size_t> i1) const;

        size_t MapToI13(std::vector<size_t> i1) const;

        std::vector<size_t> MapToI14(std::vector<size_t> i1) const;

        size_t MapToI21(std::vector<std::vector<size_t>> i2) const;
        size_t MapToI22(std::vector<std::vector<size_t>> i2) const;


        std::vector<size_t> MapToP1(std::vector<size_t> i14) const;

        size_t MapToP2() const; // Not needed for NOW!

        size_t Cfn(size_t x, size_t y) const;

        

        std::vector<std::vector<std::complex<double>>> CreateUniqueBfvs(size_t i11,
                                                                        size_t i12,
                                                                        size_t i13,
                                                                        std::vector<size_t> i14) const;

        /**
         * @brief Create the vector v_{l,m} as given in 3GPP Ts 38.214, Sec 5.2.2.2.3
         * @param l parameter l
         * @param m parameter m
         * @return the vector v_{l,m}
         */
        std::vector<std::complex<double>> CreateVecV(size_t l, size_t m) const;

        /**
         * @brief Create the vector u_m as given in 3GPP TS 38.214, Sec 5.2.2.2.3
         * @param m parameter m
         * @return the vector
         */
        std::vector<std::complex<double>> CreateVecU(size_t m) const;

        /**
         * @brief Create the Kronecker product of two vectors
         * @param vecA the first vector
         * @param vecB the second vector
         * @return a vector created by multiplying each element of vecA with the vector vecB.
         */
        static std::vector<std::complex<double>> Kroneckerproduct(
            std::vector<std::complex<double>> vecA,
            std::vector<std::complex<double>> vecB
        );

        // Derived attributes
        size_t m_o1{};  ///< Oversampling in n1-direction (typically the horizontal direction)
        size_t m_o2{};  ///< Oversampling in n2-direction (typically the vertical direction)

        size_t m_numI11{}; ///< Number of i11 values
        size_t m_numI12{}; ///< Number of i12 values
        size_t m_numI13{}; ///< Number of i13 values
        size_t m_numI14{}; ///< Number of i14 values

        std::vector<size_t> m_uniqueBfvInds;    ///< For each column in W, the beamforming vector index
        std::vector<double> m_signPhi;          ///< For each column in W, the sign before phi_n
        friend class NrPmSearchFull; 
};


} // namespace ns3

#endif // NR_CB_TYPE_TWO_H
