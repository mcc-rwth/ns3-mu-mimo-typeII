// Copyright (c) 2018 Centre Tecnologic de Telecomunicacions de Catalunya (CTTC)
//
// SPDX-License-Identifier: GPL-2.0-only

#include "nr-error-model.h"

#include "ns3/log.h"

#include "fstream"
#include "numeric"

#include "fstream"
#include "numeric"

#include <cstdlib>      // std::getenv
#include <filesystem>   // std::filesystem::create_directories
#include <fstream>
#include <map>
#include <string>
#include "ns3/simulator.h"
#include "ns3/rng-seed-manager.h"
#include "ns3/nr-eesm-error-model.h"
#include <limits>

namespace ns3
{

NS_LOG_COMPONENT_DEFINE("NrErrorModel");
NS_OBJECT_ENSURE_REGISTERED(NrErrorModel);

NrErrorModel::NrErrorModel()
    : Object()
{
    NS_LOG_FUNCTION(this);
}

NrErrorModel::~NrErrorModel()
{
}

TypeId
NrErrorModel::GetTypeId()
{
    static TypeId tid = TypeId("ns3::NrErrorModel").SetParent<Object>();
    return tid;
}

Ptr<NrErrorModelOutput>
NrErrorModel::GetTbDecodificationStatsMimo(const std::vector<MimoSinrChunk>& sinrChunks,
                                           const std::vector<int>& map,
                                           uint32_t size,
                                           uint8_t mcs,
                                           uint8_t rank,
                                           const NrErrorModelHistory& history)
{
    NS_ASSERT_MSG(!sinrChunks.empty(), "At least one SINR value is required");

    // Compute time-domain average of the SINR matrix
    auto avgSinrMat = ComputeAvgSinrMimo(sinrChunks);
    NS_ASSERT(avgSinrMat.GetNumRows() == rank);

    // Vectorize SINR matrix and convert to SpectrumValue
    /// Linearize a 2D matrix into a vector, and convert that vector to a SpectrumValue
    /// Matches layer-to-codeword mapping in TR 38.211, Table 7.3.1.3-1
    /// avgSinrMat(NrSinrMatrix) is A 2D matrix of average SINR values, dimensions nMimoLayers x
    /// nRbs vectorizedSinr(SpectrumValue) contains the (nRB * nMimoLayers) SINR values
    auto vectorizedSinr = avgSinrMat.GetVectorizedSpecVal();
    // Create a new RB map that fits the vectorized SINR values
    auto vectorizedMap = CreateVectorizedRbMap(map, rank);

    auto output = GetTbDecodificationStats(vectorizedSinr, vectorizedMap, size, mcs, history);

    if (sinrChunks[0].rnti != 0)
    {

        double sinrSum = 0.0, sinrMin = std::numeric_limits<double>::max();
        uint32_t n = 0;
        for (auto it = vectorizedSinr.ConstValuesBegin(); it != vectorizedSinr.ConstValuesEnd(); ++it)
        {
            sinrSum += *it;
            if (*it < sinrMin) sinrMin = *it;
            ++n;
        }
        double avgSinrDb = 10.0 * std::log10(sinrSum / n);
        double minSinrDb     = 10.0 * std::log10(sinrMin);   // worst RB×layer

        // EESM-effective SINR actually used by the error model (linear -> dB)
        double effSinrDb = std::numeric_limits<double>::quiet_NaN();
        if (auto eesmOut = DynamicCast<NrEesmErrorModelOutput>(output))
        {
            effSinrDb = 10.0 * std::log10(eesmOut->m_sinrEff);
        }

        const double  tNow = Simulator::Now().GetSeconds();

        // injected from main() via the --resultsDir / --runNumber / --Users args
        const char* envDir   = std::getenv("SINR_RESULTS_DIR"); // results/volume/user-8/
        const char* envRun   = std::getenv("SINR_RUN");         // 5
        const char* envUsers = std::getenv("SINR_NUSERS");      // 8
        const std::string base   = envDir   ? envDir   : "results";
        const std::string run    = envRun   ? std::string(envRun)
                                            : std::to_string(RngSeedManager::GetRun());
        const std::string nUsers = envUsers ? std::string(envUsers) : "0";

        // results/volume/user-<U>/run-<R>/sinr.csv   (filesystem::path handles the trailing slash)
        const std::filesystem::path dir  = std::filesystem::path(base) / ("run-" + run);
        const std::filesystem::path path = dir / "sinr.csv";

        // open each run's file once, keep it open
        static std::map<std::string, std::ofstream> g_sinrFiles;
        auto fit = g_sinrFiles.find(path.string());
        if (fit == g_sinrFiles.end())
        {
            std::filesystem::create_directories(dir);
            std::ofstream ofs(path, std::ios::out | std::ios::trunc);
            ofs << "time,rnti,nUsers,run,avgSinrDb,minSinrDb,effSinrDb,size,mcs,rank\n";
            fit = g_sinrFiles.emplace(path.string(), std::move(ofs)).first;
        }

        // Write data
        fit->second << tNow << "," 
                << sinrChunks[0].rnti << ","
                << nUsers << ","
                << run << ","
                << avgSinrDb << ","
                << minSinrDb << ","
                << effSinrDb << ","
                << size << ","
                << static_cast<uint32_t>(mcs) << ","
                << static_cast<uint32_t>(rank) << "\n";
    }
    
    return GetTbDecodificationStats(vectorizedSinr, vectorizedMap, size, mcs, history);
    
}

NrSinrMatrix
NrErrorModel::ComputeAvgSinrMimo(const std::vector<MimoSinrChunk>& sinrChunks)
{
    // NS_LOG_UNCOND("The Code comes here!!!!!!!!!!!!!!!!!!");
    NS_ASSERT(!sinrChunks.empty());
    if (sinrChunks.size() == 1)
    {
        // NS_LOG_UNCOND("HERE IS THE END!");
        return sinrChunks[0].mimoSinr; // Single value, no need to compute an average
    }

    // NS_LOG_UNCOND("SINR chunk size = " << sinrChunks.size());
    auto nRbs = sinrChunks[0].mimoSinr.GetNumCols();
    auto rank = sinrChunks[0].mimoSinr.GetNumRows();
    // NS_LOG_UNCOND("rank is = " << rank);
    // NS_LOG_UNCOND("Number of RBs = " << nRbs);  
    auto totDur = double{0.0};
    auto avgSinrMat = DoubleMatrixArray{rank, nRbs};
    for (const auto& chunk : sinrChunks)
    {
        const auto& sinrMat = chunk.mimoSinr;
        NS_ASSERT(sinrMat.GetNumRows() == avgSinrMat.GetNumRows());
        NS_ASSERT(sinrMat.GetNumCols() == avgSinrMat.GetNumCols());
        avgSinrMat += sinrMat * chunk.dur.GetDouble();
        totDur += chunk.dur.GetDouble();
    }

    // for (auto rb = 0; rb < nRbs; ++rb)
    // {
    //     std::cout << "SINR in the RB is a = " << avgSinrMat[0, rb] << std::endl;
    // }
    return NrSinrMatrix{avgSinrMat * (1.0 / totDur)};
}

std::vector<int>
NrErrorModel::CreateVectorizedRbMap(std::vector<int> map, uint8_t rank)
{
    auto vectorizedMap = std::vector<int>{};
    for (int iRb : map)
    {
        for (size_t layer = 0; layer < rank; layer++)
        {
            vectorizedMap.emplace_back(iRb * rank + layer);
        }
    }
    return vectorizedMap;
}

} // namespace ns3
