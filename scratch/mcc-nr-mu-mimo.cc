// Copyright (c) 2023 Centre Tecnologic de Telecomunicacions de Catalunya (CTTC)
//
// SPDX-License-Identifier: GPL-2.0-only

/**
 * @ingroup examples
 * @file cttc-nr-mimo-demo.cc
 * @brief An example that shows how to setup and use MIMO
 *
 * This example describes how to setup a simulation using MIMO. The scenario
 * consists of a simple topology, in which there
 * is only one gNB and one UE. An additional pair of gNB and UE can be enabled
 * to simulate the interference (see enableInterfNode).
 * Example creates one DL flow that goes through only BWP.
 *
 * The example prints on-screen and into the file the end-to-end result of the flow.
 * To see all the input parameters run:
 *
 * \code{.unparsed}
$ ./ns3 run cttc-nr-mimo-demo -- --PrintHelp
    \endcode
 *
 *
 * MIMO is enabled by default. To disable it run:
 *  * \code{.unparsed}
$  ./ns3 run cttc-nr-mimo-demo -- --enableMimoFeedback=0
    \endcode
 *
 *
 */


#include "ns3/antenna-module.h"
#include "ns3/applications-module.h"
#include "ns3/basic-data-calculators.h"
#include "ns3/config-store-module.h"
#include "ns3/core-module.h"
#include "ns3/fast-fading-constant-position-mobility-model.h"
#include "ns3/flow-monitor-module.h"
#include "ns3/internet-apps-module.h"
#include "ns3/internet-module.h"
#include "ns3/mobility-module.h"
#include "ns3/network-module.h"
#include "ns3/nr-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/stats-module.h"
#include "ns3/traffic-generator-helper.h"
#include "ns3/udp-server.h"

#include <map>
#include <fstream>
#include <vector>
#include <filesystem>

using namespace ns3;
NS_LOG_COMPONENT_DEFINE("MccNrMuMimoDemo");


struct UeAppTraceInfo
{
    uint32_t ueIdx;
    uint16_t rnti;
    Ipv4Address ip;
    Ptr<UdpServer> server;
    uint64_t lastRxPackets;
};

// std::map<uint32_t, uint64_t> g_rxBytes;
static std::ofstream g_thrFile;
uint32_t g_runNumber = 1;
std::string g_resultsDir = "out";
static std::vector<UeAppTraceInfo> g_ueTraceInfo;
static uint32_t g_udpPacketSize = 0;

static std::map<uint64_t, std::unique_ptr<std::ofstream>> g_dlSinrFiles;

static void
DlDataSinrTracePerUe(uint32_t ueIdx,
                     uint64_t imsi,
                     uint16_t cellId,
                     uint16_t rnti,
                     double sinrLin,
                     uint16_t bwpId)
{
    std::string runDir = g_resultsDir + "/run-" + std::to_string(g_runNumber);
    std::filesystem::create_directories(runDir);

    auto& fp = g_dlSinrFiles[imsi];
    if (!fp)
    {
        std::string filename = runDir + "/dl_sinr_ue" + std::to_string(ueIdx) +
                               "_imsi" + std::to_string(imsi) + ".csv";
        fp = std::make_unique<std::ofstream>(filename, std::ios::out | std::ios::trunc);

        (*fp) << "time_s,run,ueIdx,imsi,cellId,rnti,bwpId,sinrLin,sinrDb\n";
        fp->setf(std::ios::fixed);
        fp->precision(9);
    }

    double sinrDb = 10.0 * std::log10(std::max(sinrLin, 1e-15));

    (*fp) << Simulator::Now().GetSeconds() << ","
          << g_runNumber << ","
          << ueIdx << ","
          << imsi << ","
          << cellId << ","
          << rnti << ","
          << bwpId << ","
          << sinrLin << ","
          << sinrDb << "\n";
}

static void
ConnectDlSinrTraces(const NetDeviceContainer& ueNetDev)
{
    for (uint32_t i = 0; i < ueNetDev.GetN(); ++i)
    {
        Ptr<NrUeNetDevice> ue = DynamicCast<NrUeNetDevice>(ueNetDev.Get(i));
        NS_ABORT_MSG_IF(!ue, "Device is not NrUeNetDevice");

        uint64_t imsi = ue->GetImsi();

        // If you use one BWP, GetPhy(0) is enough.
        // If you use multiple BWPs, loop over them.
        Ptr<NrUePhy> uePhy = ue->GetPhy(0);
        NS_ABORT_MSG_IF(!uePhy, "UE PHY is null");

        uePhy->TraceConnectWithoutContext(
            "DlDataSinr",
            MakeBoundCallback(&DlDataSinrTracePerUe, i, imsi));
    }
}

static void
CloseDlSinrFiles()
{
    for (auto& [imsi, fp] : g_dlSinrFiles)
    {
        if (fp && fp->is_open())
        {
            fp->close();
        }
    }
    g_dlSinrFiles.clear();
}

void
LogPerUeThroughput()
{
    double now = ns3::Simulator::Now().GetSeconds();
    double windowSec = 0.01;

    for (auto& ue : g_ueTraceInfo)
    {
        uint64_t curRxPackets = ue.server->GetReceived();
        uint64_t diffPackets = curRxPackets - ue.lastRxPackets;
        ue.lastRxPackets = curRxPackets;

        uint64_t diffBytes = diffPackets * g_udpPacketSize;
        double thrMbps = (diffBytes * 8.0) / 1e6 / windowSec;

        g_thrFile << g_runNumber << ","
                  << now << ","
                  << ue.ueIdx << ","
                  << ue.rnti << ","
                  << ue.ip << ","
                  << thrMbps << ","
                  << curRxPackets << ","
                  << (curRxPackets * g_udpPacketSize)
                  << "\n";
    }

    ns3::Simulator::Schedule(ns3::MilliSeconds(10), &LogPerUeThroughput);

}



struct CqiFeedbackTraceStats
{
    Ptr<MinMaxAvgTotalCalculator<uint8_t>> m_ri;
    Ptr<MinMaxAvgTotalCalculator<uint8_t>> m_mcs;

    CqiFeedbackTraceStats()
    {
        m_ri = CreateObject<MinMaxAvgTotalCalculator<uint8_t>>();
        m_mcs = CreateObject<MinMaxAvgTotalCalculator<uint8_t>>();
    }

    CqiFeedbackTraceStats(uint8_t rank, uint8_t mcs)
    {
        m_ri = CreateObject<MinMaxAvgTotalCalculator<uint8_t>>();
        m_ri->Update(rank);
        m_mcs = CreateObject<MinMaxAvgTotalCalculator<uint8_t>>();
        m_mcs->Update(mcs);
    }
};

void
CqiFeedbackTracedCallback(std::map<uint16_t, CqiFeedbackTraceStats>* stats,
                          uint16_t rnti,
                          [[maybe_unused]] uint8_t cqi,
                          uint8_t mcs,
                          uint8_t rank)
{
    auto it = stats->find(rnti);
    if (it != stats->end())
    {
        it->second.m_ri->Update(rank);
        it->second.m_mcs->Update(mcs);
    }
    else
    {
        (*stats)[rnti] = CqiFeedbackTraceStats(rank, mcs);
    }
}

void
LogUePositions(const NodeContainer& ueNodes,
               std::string filename,
               Time interval)
{
    std::ofstream file(filename, std::ios::app);

    if (!file.is_open())
    {
        NS_LOG_UNCOND("Error opening file " << filename);
        return;
    }

    double now = Simulator::Now().GetSeconds();

    for (int i = 0; i < ueNodes.GetN(); ++i)
    {
        Ptr<Node> node = ueNodes.Get(i);
        Ptr<MobilityModel> mob = node->GetObject<MobilityModel>();

        if (!mob)
        {
            continue;
        }

        Vector pos = mob->GetPosition();

        file << now << ","        // time
             << i << ","          // UE index
             << pos.x << ","
             << pos.y << ","
             << pos.z << "\n";
    }

    file.close();

    // Schedule next logging
    Simulator::Schedule(interval, &LogUePositions, ueNodes, filename, interval);
}

int
main(int argc, char* argv[])
{
    
    uint32_t runNumber = 1;
    auto startExecTime = std::chrono::system_clock::now();
    bool enableMimoFeedback = true;
    bool useConfigSetDefault = false;
    uint8_t csiFlags = 3;

    NrHelper::AntennaParams apUe;
    NrHelper::AntennaParams apGnb;
    apUe.antennaElem = "ns3::ThreeGppAntennaModel";
    apUe.nAntCols = 2;
    apUe.nAntRows = 2;
    apUe.nHorizPorts = 2;
    apUe.nVertPorts = 1;
    apUe.isDualPolarized = false;
    apGnb.antennaElem = "ns3::ThreeGppAntennaModel";
    apGnb.nAntCols = 4;
    apGnb.nAntRows = 4;
    apGnb.nHorizPorts = 4;
    apGnb.nVertPorts = 4;
    apGnb.isDualPolarized = false;
    double downtiltAngleGnb = 10;

    // The polarization slant angle in degrees in case of x-polarized
    double polSlantAngleGnb = 0.0;
    double polSlantAngleUe = 90.0;
    // The bearing angles in degrees
    double bearingAngleGnb = 0.0;
    double bearingAngleUe = 180.0;

    std::string trafficType = "cbr";
    // Traffic parameters
    uint32_t udpPacketSize = 1000;
    // Packet interval is 40000 ns to reach 200 Mbps
    // For MCS Table 2, and 10 MHz BW, 200 Mbps can be achieved by using 4 MIMO streams
    Time packetInterval = MilliSeconds(30);
    Time udpAppStartTime = MilliSeconds(200);

    // Interference
    bool enableInterfNode = false; // if true an additional pair of gNB and UE will be created to
                                   // create an interference towards the original pair
    double interfDistance =
        1000.0; // the distance in meters between the gNB1 and the interfering gNB2
    double interfPolSlantDelta = 0; // the difference between the pol. slant angle between the
                                    // original node and the interfering one

    // Other simulation scenario parameters
    Time simTime = MilliSeconds(10000);
    uint16_t gnbUeDistance = 20; // meters
    uint16_t numerology = 0;
    double centralFrequency = 3.5e9;
    double bandwidth = 10e6;
    double txPowerGnb = 28; // dBm
    double txPowerUe = 23;  // dBm
    uint16_t updatePeriodMs = 50;
    std::string errorModel = "ns3::NrEesmIrT2";
    std::string scheduler = "ns3::NrMacSchedulerTdmaRR";
    std::string beamformingMethod = "ns3::DirectPathBeamforming";

    uint32_t wbPmiUpdateIntervalMs = 1; // Wideband PMI update interval in ms
    uint32_t sbPmiUpdateIntervalMs = 1;  // Subband PMI update interval in ms

    // Default channel condition
    std::string losCondition = "Default";
    NrHelper::MimoPmiParams mimoPmiParams;
    mimoPmiParams.subbandSize = 8;
    double xyVelocity = 0;

    // Number of users in the simulation
    uint16_t ueNumber = 2;
    double sdmaThres = 3e-1;
    uint32_t maxSdmaUsers = 4;
    bool csiBaseline = true;
    uint32_t maxStreams = 4;

    // Where the example stores the output files.
    std::string simTag = "default";
    std::string outputDir = "./out";
    bool logging = false;
    std::string resultsDir = "out";

    CommandLine cmd(__FILE__);
    /**
     * The main parameters for testing MIMO
     */
    cmd.AddValue("enableMimoFeedback", "Enables MIMO feedback", enableMimoFeedback);
    cmd.AddValue(
        "pmSearchMethod",
        "Precoding matrix search method, currently implemented only exhaustive search method"
        "[ns3::NrPmSearchFull, ns3::NrPmSearchFast, ns3::NrPmSearchIdeal, ns3::NrPmSearchSasaoka, "
        "ns3::NrPmSearchMaleki (requires extra dependencies)]",
        mimoPmiParams.pmSearchMethod);
    cmd.AddValue("fullSearchCb",
                 "The codebook to be used for the full search. Available codebooks are "
                 "a) ns3::NrCbTwoPort, the two-port codebook defined in 3GPP TS 38.214 Table "
                 "5.2.2.2.1-1, and"
                 "b) ns3::NrCbTypeOneSp, Type-I Single-Panel Codebook 3GPP TS 38.214 Rel. 15, "
                 "Sec. 5.2.2.2.1 supporting codebook mode 1 only, and limited to rank 4.",
                 mimoPmiParams.fullSearchCb);
    cmd.AddValue("rankLimit", "The maximum rank number to be used.", mimoPmiParams.rankLimit);
    cmd.AddValue("rankTechnique",
                 "Technique used for RI selection by Fast and Sasaoka PMI selection [SVD, "
                 "WaterFilling, Sasaoka]",
                 mimoPmiParams.rankTechnique);
    cmd.AddValue("rankThreshold", "Threshold used by rankTechnique", mimoPmiParams.rankThreshold);
    cmd.AddValue("subbandSize", "Sub-band size for downsampling", mimoPmiParams.subbandSize);
    cmd.AddValue("downsamplingTechnique",
                 "Sub-band downsampling technique",
                 mimoPmiParams.downsamplingTechnique);
    cmd.AddValue("numRowsGnb", "Number of antenna rows at the gNB", apGnb.nAntRows);
    cmd.AddValue("numRowsUe", "Number of antenna rows at the UE", apUe.nAntRows);
    cmd.AddValue("numColumnsGnb", "Number of antenna columns at the gNB", apGnb.nAntCols);
    cmd.AddValue("numColumnsUe", "Number of antenna columns at the UE", apUe.nAntCols);
    cmd.AddValue("numVPortsGnb",
                 "Number of vertical ports of the antenna at the gNB",
                 apGnb.nVertPorts);
    cmd.AddValue("numVPortsUe",
                 "Number of vertical ports of the antenna at the UE",
                 apUe.nVertPorts);
    cmd.AddValue("numHPortsGnb",
                 "Number of horizontal ports of the antenna the gNB",
                 apGnb.nHorizPorts);
    cmd.AddValue("numHPortsUe",
                 "Number of horizontal ports of the antenna at the UE",
                 apUe.nHorizPorts);
    cmd.AddValue("xPolGnb",
                 "Whether the gNB antenna array has the cross polarized antenna "
                 "elements.",
                 apGnb.isDualPolarized);
    cmd.AddValue("xPolUe",
                 "Whether the UE antenna array has the cross polarized antenna "
                 "elements.",
                 apUe.isDualPolarized);
    cmd.AddValue("polSlantAngleGnb",
                 "Polarization slant angle of gNB in degrees",
                 polSlantAngleGnb);
    cmd.AddValue("polSlantAngleUe", "Polarization slant angle of UE in degrees", polSlantAngleUe);
    cmd.AddValue("bearingAngleGnb", "Bearing angle of gNB in degrees", bearingAngleGnb);
    cmd.AddValue("bearingAngleUe", "Bearing angle of UE in degrees", bearingAngleUe);
    cmd.AddValue("downtiltAngleGnb", "Downtilt angle of gNB in degrees", downtiltAngleGnb);
    cmd.AddValue("enableInterfNode", "Whether to enable an interfering node", enableInterfNode);
    cmd.AddValue("wbPmiUpdateInterval",
                 "Wideband PMI update interval in ms",
                 wbPmiUpdateIntervalMs);
    cmd.AddValue("sbPmiUpdateInterval", "Subband PMI update interval in ms", sbPmiUpdateIntervalMs);
    cmd.AddValue("interfDistance",
                 "The distance between the gNB1 and the interfering gNB2 (the original and the "
                 "interfering one)",
                 interfDistance);
    cmd.AddValue("interfPolSlantDelta",
                 "The difference between the pol. slant angles of the original pairs of gNB and UE "
                 "and the interfering one",
                 interfPolSlantDelta);
    cmd.AddValue("csiFlags", "CsiFlags to be configured. See NrHelper::CsiFlags", csiFlags);
    /**
     * Other simulation parameters
     */
    cmd.AddValue("trafficType",
                 "Traffic type to be installed at the source: cbr or ftp.",
                 trafficType);
    cmd.AddValue("packetSize",
                 "packet size in bytes to be used by best effort traffic",
                 udpPacketSize);
    cmd.AddValue("packetInterval", "Inter-packet interval for CBR traffic", packetInterval);
    cmd.AddValue("simTime", "Simulation time", simTime);
    cmd.AddValue("numerology", "The numerology to be used", numerology);
    cmd.AddValue("centralFrequency", "The system frequency to be used in band 1", centralFrequency);
    cmd.AddValue("bandwidth", "The system bandwidth to be used", bandwidth);
    cmd.AddValue("txPowerGnb", "gNB TX power", txPowerGnb);
    cmd.AddValue("txPowerUe", "UE TX power", txPowerUe);
    cmd.AddValue("gnbUeDistance",
                 "The distance between the gNB and the UE in the scenario",
                 gnbUeDistance);
    cmd.AddValue(
        "updatePeriodMs",
        "Channel update period in ms. If set to 0 then the channel update will be disabled",
        updatePeriodMs);
    cmd.AddValue("errorModel",
                 "Error model: ns3::NrEesmCcT1, ns3::NrEesmCcT2, "
                 "ns3::NrEesmIrT1, ns3::NrEesmIrT2, ns3::NrLteMiErrorModel",
                 errorModel);
    cmd.AddValue("scheduler",
                 "The scheduler: ns3::NrMacSchedulerTdmaRR, "
                 "ns3::NrMacSchedulerTdmaPF, ns3::NrMacSchedulerTdmaMR,"
                 "ns3::NrMacSchedulerTdmaQos, ns3::NrMacSchedulerOfdmaRR, "
                 "ns3::NrMacSchedulerOfdmaPF, ns3::NrMacSchedulerOfdmaMR,"
                 "ns3::NrMacSchedulerOfdmaQos",
                 scheduler);
    cmd.AddValue("beamformingMethod",
                 "The beamforming method: ns3::CellScanBeamforming,"
                 "ns3::CellScanQuasiOmniBeamforming,"
                 "ns3::DirectPathBeamforming,"
                 "ns3::QuasiOmniDirectPathBeamforming,"
                 "ns3::DirectPathQuasiOmniBeamforming,"
                 "ns3::KronBeamforming,"
                 "ns3::KronQuasiOmniBeamforming",
                 beamformingMethod);
    cmd.AddValue("losCondition",
                 "Default - for 3GPP channel condition model,"
                 "LOS - for always LOS channel condition model,"
                 "NLOS - for always NLOS channel condition model",
                 losCondition);
    cmd.AddValue("simTag",
                 "tag to be appended to output filenames to distinguish simulation campaigns",
                 simTag);
    cmd.AddValue("outputDir", "directory where to store simulation results", outputDir);
    cmd.AddValue("logging", "Enable logging", logging);
    cmd.AddValue("useConfigSetDefault",
                 "Configure via Config::SetDefault instead of the MimoPmiParams structure",
                 useConfigSetDefault);
    cmd.AddValue("xyVelocity",
                 "Velocity in X and Y directions m/s for fake fading model.",
                 xyVelocity);
    cmd.AddValue("Users",
                 "The number of the users in the simulation.",
                 ueNumber);
    cmd.AddValue("threshold",
                 "Threshold of SDMA grouping",
                 sdmaThres);
    cmd.AddValue("maxSdmaUsers",
                    "Maximum number of users in a group",
                    maxSdmaUsers);
    cmd.AddValue("FullCsiBaseline",
                    "Baseline SDMA with full CSI",
                    csiBaseline);
    cmd.AddValue("MaxStreams",
                    "Maximum number of streams",
                    maxStreams);
    cmd.AddValue("runNumber",
                 "Independent run number",
                 runNumber);
    cmd.AddValue("resultsDir",
                 "Directory where results are stored",
                 resultsDir);
    // Parse the command line
    cmd.Parse(argc, argv);

    setenv("SINR_RESULTS_DIR", resultsDir.c_str(), 1);
    setenv("SET_RUN", std::to_string(runNumber).c_str(), 1);
    setenv("SINR_NUSERS", std::to_string(ueNumber).c_str(), 1);

    RngSeedManager::SetSeed(42);
    RngSeedManager::SetRun(runNumber);

    // convert angle values into radians
    apUe.bearingAngle = bearingAngleUe * (M_PI / 180);
    apUe.polSlantAngle = polSlantAngleUe * (M_PI / 180);
    apGnb.bearingAngle = bearingAngleGnb * (M_PI / 180);
    apGnb.polSlantAngle = polSlantAngleGnb * (M_PI / 180);

    NS_ABORT_IF(centralFrequency < 0.5e9 && centralFrequency > 100e9);

    if (logging)
    {
        LogComponentEnable("UdpClient", LOG_LEVEL_INFO);
        LogComponentEnable("UdpServer", LOG_LEVEL_INFO);
        LogComponentEnable("NrPdcp", LOG_LEVEL_INFO);
    }

    Config::SetDefault("ns3::NrRlcUm::MaxTxBufferSize", UintegerValue(999999999));
    Config::SetDefault("ns3::ThreeGppChannelModel::UpdatePeriod",
                       TimeValue(MilliSeconds(updatePeriodMs)));

    uint16_t pairsToCreate = 1;
    

    NodeContainer gnbContainer;
    gnbContainer.Create(pairsToCreate);
    NodeContainer ueContainer;
    ueContainer.Create(ueNumber);

    /**
     * We configure the mobility model to ConstantPositionMobilityModel.
     * The default topology is the following:
     *
     * gNB1.................UE1................UE2........................gNB2(interferer)
     *(0.0, 0.0, 25.0)  (d, 0.0, 1.5)    (interfDistance/2, 0.0, 1.5)    (interfDistance,0.0, 25.0)
     * bearingAngle=0   bearingAngle=180 bearingAngle=0                   bearingAngle=180
     */
    MobilityHelper gnbMobility;
    gnbMobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    Ptr<ListPositionAllocator> gnbPositionAlloc = CreateObject<ListPositionAllocator>();
    gnbPositionAlloc->Add(Vector(0.0, 0.0, 10.0));

    // MobilityHelper ueMobilityStatic;

    // ueMobilityStatic.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    // Ptr<ListPositionAllocator> uePositionAlloc1 = CreateObject<ListPositionAllocator>();

    // for (uint8_t ue = 0; ue < ueNumber; ++ue)
    // {
    //     uePositionAlloc1->Add(Vector((std::pow(1, ue) * (ue+1) * gnbUeDistance), 0.0, 1.5));
    // }

    MobilityHelper ueMobility;
    uint16_t distMin = 50, distMax = 250;

    // Optional: initial positions (can also skip this)
    ueMobility.SetPositionAllocator("ns3::RandomRectanglePositionAllocator",
                                "X", StringValue("ns3::UniformRandomVariable[Min="+std::to_string(distMin)+"|Max="+std::to_string(distMax)+"]"),
                                "Y", StringValue("ns3::UniformRandomVariable[Min="+std::to_string(distMin)+"|Max="+std::to_string(distMax)+"]"),
                                "Z", DoubleValue(1.5));

    // Random walk within bounded area
    // ueMobility.SetMobilityModel("ns3::RandomWalk2dMobilityModel",
    //                         "Bounds", RectangleValue(Rectangle(distMin, distMax, distMin, distMax)),
    //                         "Speed", StringValue("ns3::ConstantRandomVariable[Constant=0.0]"),
    //                         "Distance", DoubleValue(0.0));



    gnbMobility.SetPositionAllocator(gnbPositionAlloc);
    // ueMobilityStatic.SetPositionAllocator(uePositionAlloc1);
    // ueMobilityMoving.SetPositionAllocator(uePositionAlloc2);

    gnbMobility.Install(gnbContainer.Get(0));
    ueMobility.Install(ueContainer);
    // ueMobilityStatic.Install(ueContainer);


    for (int ue = 0; ue <ueContainer.GetN(); ++ue)
    {
        Ptr<Node> ueNode = ueContainer.Get(ue);
        Ptr<MobilityModel> mob = ueNode->GetObject<MobilityModel>();

        NS_LOG_UNCOND(mob->GetPosition().x << ", " << mob->GetPosition().y << ", " << mob->GetPosition().z);
    }

    std::string filenameUe = resultsDir + "/run-" + std::to_string(runNumber) + "/ue-positions.csv";

    // Optional: write header once
    std::ofstream file(filenameUe);
    file << "time,ue,x,y,z\n";
    file.close();

    // Start periodic logging (e.g., every 0.01s = 10 ms)
    Time logInterval = MilliSeconds(100);

    Simulator::Schedule(logInterval,
                        &LogUePositions,
                        ueContainer,
                        filenameUe,
                        logInterval);


    /**
     * Create the NR helpers that will be used to create and setup NR devices, spectrum, ...
     */
    Ptr<NrPointToPointEpcHelper> nrEpcHelper = CreateObject<NrPointToPointEpcHelper>();
    Ptr<IdealBeamformingHelper> idealBeamformingHelper = CreateObject<IdealBeamformingHelper>();
    Ptr<NrHelper> nrHelper = CreateObject<NrHelper>();
    nrHelper->SetBeamformingHelper(idealBeamformingHelper);
    nrHelper->SetEpcHelper(nrEpcHelper);
    /**
     * Prepare spectrum. Prepare one operational band, containing
     * one component carrier, and a single bandwidth part
     * centered at the frequency specified by the input parameters.
     *
     *
     * The configured spectrum division is:
     * ------------Band--------------
     * ------------CC1----------------
     * ------------BWP1---------------
     */

    CcBwpCreator ccBwpCreator;
    const uint8_t numCcPerBand = 1;
    CcBwpCreator::SimpleOperationBandConf bandConf(centralFrequency, bandwidth, numCcPerBand);
    OperationBandInfo band = ccBwpCreator.CreateOperationBandContiguousCc(bandConf);
    // Create the channel helper
    Ptr<NrChannelHelper> channelHelper = CreateObject<NrChannelHelper>();
    // Set the channel using the scenario and user input
    channelHelper->ConfigureFactories("UMi", losCondition, "ThreeGpp");
    // Set the channel update period and shadowing
    channelHelper->SetChannelConditionModelAttribute("UpdatePeriod",
                                                     TimeValue(MilliSeconds(updatePeriodMs)));
    channelHelper->SetPathlossAttribute("ShadowingEnabled", BooleanValue(true));
    // Create and set the channel with the band
    channelHelper->AssignChannelsToBands({band});

    // Configure NrHelper, prepare most of the parameters that will be used in the simulation.
    nrHelper->SetAttribute("CsiFeedbackFlags", UintegerValue(csiFlags));
    nrHelper->SetDlErrorModel(errorModel);
    nrHelper->SetUlErrorModel(errorModel);
    nrHelper->SetGnbDlAmcAttribute("AmcModel", EnumValue(NrAmc::ErrorModel));
    nrHelper->SetGnbUlAmcAttribute("AmcModel", EnumValue(NrAmc::ErrorModel));
    nrHelper->SetSchedulerTypeId(TypeId::LookupByName(scheduler));
    idealBeamformingHelper->SetAttribute("BeamformingMethod",
                                         TypeIdValue(TypeId::LookupByName(beamformingMethod)));
    nrHelper->SetSchedulerAttribute("compatThres", DoubleValue(sdmaThres));
    nrHelper->SetSchedulerAttribute("maxSdmaUsers", UintegerValue(maxSdmaUsers));
    nrHelper->SetSchedulerAttribute("FullCsiBaseline", BooleanValue(csiBaseline));
    nrHelper->SetSchedulerAttribute("MaxStreams", UintegerValue(maxStreams));
    // Core latency
    nrEpcHelper->SetAttribute("S1uLinkDelay", TimeValue(MilliSeconds(0)));

    // We can configure via Config::SetDefault
    if (enableMimoFeedback)
    {
        // We can configure not only via Config::SetDefault, but also via the MimoPmiParams
        // structure
        if (useConfigSetDefault)
        {
            Config::SetDefault("ns3::NrHelper::EnableMimoFeedback", BooleanValue(true));
            Config::SetDefault("ns3::NrPmSearch::SubbandSize", UintegerValue(16));
        }
        else
        {
            nrHelper->SetupMimoPmi(mimoPmiParams);
        }
    }

    /**
     * Configure gNb antenna
     */
    nrHelper->SetupGnbAntennas(apGnb);
    // TODO consider adding DowntiltAngle to AntennaParams
    nrHelper->SetGnbAntennaAttribute("DowntiltAngle", DoubleValue(downtiltAngleGnb * M_PI / 180.0));
    /**
     * Configure UE antenna
     */
    nrHelper->SetupUeAntennas(apUe);

    nrHelper->SetGnbPhyAttribute("Numerology", UintegerValue(numerology));
    nrHelper->SetGnbPhyAttribute("TxPower", DoubleValue(txPowerGnb));
    nrHelper->SetUePhyAttribute("TxPower", DoubleValue(txPowerUe));
    nrHelper->SetUePhyAttribute("WbPmiUpdateInterval",
                                TimeValue(MilliSeconds(wbPmiUpdateIntervalMs)));
    nrHelper->SetUePhyAttribute("SbPmiUpdateInterval",
                                TimeValue(MilliSeconds(sbPmiUpdateIntervalMs)));

    

    uint32_t bwpId = 0;
    // gNb routing between bearer type and bandwidth part
    nrHelper->SetGnbBwpManagerAlgorithmAttribute("GBR_CONV_VOICE", UintegerValue(bwpId));
    // UE routing between bearer type and bandwidth part
    nrHelper->SetUeBwpManagerAlgorithmAttribute("GBR_CONV_VOICE", UintegerValue(bwpId));
    BandwidthPartInfoPtrVector allBwps;
    allBwps = CcBwpCreator::GetAllBwps({band});

    /**
     * Finally, create the gNB and the UE device.
     */
    NetDeviceContainer gnbNetDev = nrHelper->InstallGnbDevice(gnbContainer, allBwps);
    NetDeviceContainer ueNetDev = nrHelper->InstallUeDevice(ueContainer, allBwps);

    if (enableInterfNode)
    {
        nrHelper->GetGnbPhy(gnbNetDev.Get(1), 0)
            ->GetSpectrumPhy()
            ->GetAntenna()
            ->SetAttribute("BearingAngle", DoubleValue(0));
        nrHelper->GetUePhy(ueNetDev.Get(1), 0)
            ->GetSpectrumPhy()
            ->GetAntenna()
            ->SetAttribute("BearingAngle", DoubleValue(M_PI));
        if (interfPolSlantDelta)
        {
            // reconfigure the polarization slant angle of the interferer
            nrHelper->GetGnbPhy(gnbNetDev.Get(1), 0)
                ->GetSpectrumPhy()
                ->GetAntenna()
                ->SetAttribute(
                    "PolSlantAngle",
                    DoubleValue((polSlantAngleGnb + interfPolSlantDelta) * (M_PI / 180)));
            nrHelper->GetUePhy(ueNetDev.Get(1), 0)
                ->GetSpectrumPhy()
                ->GetAntenna()
                ->SetAttribute("PolSlantAngle",
                               DoubleValue((polSlantAngleUe + interfPolSlantDelta) * (M_PI / 180)));
        }
    }

    /**
     * Fix the random stream throughout the nr, propagation, and spectrum
     * modules classes. This configuration is extremely important for the
     * reproducibility of the results.
     */
    int64_t randomStream = 1;
    randomStream += nrHelper->AssignStreams(gnbNetDev, randomStream);
    randomStream += nrHelper->AssignStreams(ueNetDev, randomStream);

    std::map<uint16_t, CqiFeedbackTraceStats> cqiTraces;
    for (auto it = ueNetDev.Begin(); it != ueNetDev.End(); ++it)
    {
        auto cqiCb = MakeBoundCallback(&CqiFeedbackTracedCallback, &cqiTraces);
        nrHelper->GetUePhy(*it, 0)->TraceConnectWithoutContext("CqiFeedbackTrace", cqiCb);
    }

    // create the Internet and install the IP stack on the UEs
    // get SGW/PGW and create a single RemoteHost
    auto [remoteHost, remoteHostIpv4Address] =
        nrEpcHelper->SetupRemoteHost("100Gb/s", 2500, Seconds(0.000));

    InternetStackHelper internet;
    internet.Install(ueContainer);
    Ipv4InterfaceContainer ueIpIface =
        nrEpcHelper->AssignUeIpv4Address(NetDeviceContainer(ueNetDev));

    // attach each UE to its gNB according to desired scenario
    // nrHelper->AttachToGnb(ueNetDev.Get(0), gnbNetDev.Get(0));
    // nrHelper->AttachToGnb(ueNetDev.Get(1), gnbNetDev.Get(0));
    // nrHelper->AttachToGnb(ueNetDev.Get(2), gnbNetDev.Get(0));
    for (size_t users = 0; users < ueNetDev.GetN(); ++users)
    {
        nrHelper->AttachToGnb(ueNetDev.Get(users), gnbNetDev.Get(0));
    }

    /**
     * Install DL traffic part.
     */
    uint16_t dlPort = 1234;
    ApplicationContainer serverApps;
    // The sink will always listen to the specified ports
    UdpServerHelper dlPacketSink(dlPort);
    // The server, that is the application which is listening, is installed in the UE
    serverApps.Add(dlPacketSink.Install(ueContainer));

    // The bearer that will carry the traffic
    NrEpsBearer epsBearer(NrEpsBearer::GBR_CONV_VOICE);

    // The filter for the traffic
    Ptr<NrEpcTft> dlTft = Create<NrEpcTft>();
    NrEpcTft::PacketFilter dlPktFilter;
    dlPktFilter.localPortStart = dlPort;
    dlPktFilter.localPortEnd = dlPort;
    dlTft->Add(dlPktFilter);

    /**
     * Let's install the applications!
     */

    ApplicationContainer clientApps;

    if (trafficType == "cbr")
    {
        // Option A: same dlPort for all UEs (OK if each UE has its own IP)
        // Option B (recommended): unique port per UE to make tracing/debugging easier
        uint16_t baseDlPort = dlPort;

        for (uint32_t u = 0; u < ueNetDev.GetN(); ++u)
        {
            UdpClientHelper dlClient;

            dlClient.SetAttribute("MaxPackets", UintegerValue(0xFFFFFFFF));
            dlClient.SetAttribute("PacketSize", UintegerValue(udpPacketSize));
            dlClient.SetAttribute("Interval", TimeValue(packetInterval));

            // uint16_t uePort = baseDlPort + u; // unique per UE (recommended)
            uint16_t uePort = baseDlPort; // unique per UE (recommended)

            dlClient.SetAttribute(
                "Remote",
                AddressValue(addressUtils::ConvertToSocketAddress(ueIpIface.GetAddress(u), uePort)));

            clientApps.Add(dlClient.Install(remoteHost));

            // Activate a dedicated bearer for THIS UE
            nrHelper->ActivateDedicatedEpsBearer(ueNetDev.Get(u), epsBearer, dlTft);
        }
    }

    g_ueTraceInfo.clear();
    g_udpPacketSize = udpPacketSize;

    for (uint32_t u = 0; u < ueContainer.GetN(); ++u)
    {
        Ptr<UdpServer> server = DynamicCast<UdpServer>(serverApps.Get(u));
        NS_ASSERT(server);

        Ptr<NetDevice> dev = ueNetDev.Get(u);
        Ptr<NrUeNetDevice> nrUeDev = DynamicCast<NrUeNetDevice>(dev);
        NS_ASSERT(nrUeDev);

        uint16_t rnti = nrUeDev->GetRrc()->GetRnti();
        Ipv4Address ip = ueIpIface.GetAddress(u);

        g_ueTraceInfo.push_back({u, rnti, ip, server, 0});
    }

    g_runNumber = runNumber;
    g_resultsDir = resultsDir;

    std::filesystem::create_directories(g_resultsDir);

    // g_thrFile.open("throughput-time.csv", std::ofstream::out | std::ofstream::trunc);
    std::ostringstream fname;
    fname << g_resultsDir << "/throughput-time-run-"
        << std::setw(3) << std::setfill('0') << g_runNumber
        << ".csv";

    g_thrFile.open(fname.str(),
               std::ofstream::out | std::ofstream::trunc);
    g_thrFile << "run,Time(s),ueIdx,rnti,ip,throughputMbps,totalRxPackets,TotalRxBytes\n";

    Simulator::Schedule(MilliSeconds(10), &LogPerUeThroughput);

    ConnectDlSinrTraces(ueNetDev);

    Ptr<NrPhyRxTrace> phyTrace = CreateObject<NrPhyRxTrace>();
    phyTrace->SetResultsFolder(resultsDir + "/run-" + std::to_string(runNumber) + "/");
    phyTrace->SetSimTag("_run" + std::to_string(runNumber));

    Ptr<NrMacSchedulingStats> schedTrace = CreateObject<NrMacSchedulingStats>();
    schedTrace->SetAttribute("DlOutputFilename", StringValue(resultsDir + "/run-" + std::to_string(runNumber) + "/" + "DlSchedStats.txt"));

    // DL PHY trace (RxPacketTrace, at the UE)
    Config::Connect(
        "/NodeList/*/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/NrUePhy/NrSpectrumPhy/RxPacketTraceUe",
        MakeBoundCallback(&NrPhyRxTrace::RxPacketTraceUeCallback, phyTrace));

    // DL MAC scheduling stats
    Config::Connect(
        "/NodeList/*/DeviceList/*/$ns3::NrGnbNetDevice/BandwidthPartMap/*/NrGnbMac/DlScheduling",
        MakeBoundCallback(&NrMacSchedulingStats::DlSchedulingCallback, schedTrace));

    // start UDP server and client apps
    serverApps.Start(udpAppStartTime);
    clientApps.Start(udpAppStartTime);
    serverApps.Stop(simTime);
    clientApps.Stop(simTime);

    // enable the traces provided by the nr module
    // nrHelper->EnableTraces();

    FlowMonitorHelper flowmonHelper;
    NodeContainer endpointNodes;
    endpointNodes.Add(remoteHost);
    endpointNodes.Add(ueContainer);

    Ptr<ns3::FlowMonitor> monitor = flowmonHelper.Install(endpointNodes);
    monitor->SetAttribute("DelayBinWidth", DoubleValue(0.001));
    monitor->SetAttribute("JitterBinWidth", DoubleValue(0.001));
    monitor->SetAttribute("PacketSizeBinWidth", DoubleValue(20));

    Simulator::Stop(simTime);
    Simulator::Run();

    // Print per-flow statistics
    monitor->CheckForLostPackets();
    Ptr<Ipv4FlowClassifier> classifier =
        DynamicCast<Ipv4FlowClassifier>(flowmonHelper.GetClassifier());
    FlowMonitor::FlowStatsContainer stats = monitor->GetFlowStats();

    double averageFlowThroughput = 0.0;
    double averageFlowDelay = 0.0;

    std::ofstream outFile;
    std::string filename = outputDir + "/" + simTag;
    outFile.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);
    if (!outFile.is_open())
    {
        std::cerr << "Can't open file " << filename << std::endl;
        return 1;
    }

    outFile.setf(std::ios_base::fixed);

    

    // calculate the execution time
    auto endExecTime = std::chrono::system_clock::now();
    double execTimeSec = std::chrono::duration<double>(endExecTime - startExecTime).count();

    std::ofstream execFile;
    execFile.open(resultsDir + "/run-" + std::to_string(runNumber) + "/execution-time.csv", std::ios::out);

    execFile << "RunId,ExecutionTimeSeconds\n";
    execFile << runNumber << "," << execTimeSec << "\n";

    execFile.close();


    double rank;
    double mcs;

    double averageMcsForAllUes = 0.0;
    double averageRiForAllUes = 0.0;
    for (const auto& ue : cqiTraces)
    {
        averageRiForAllUes += ue.second.m_ri->getMean();
        averageMcsForAllUes += ue.second.m_mcs->getMean();
    }

    if (ueNetDev.GetN() != cqiTraces.size())
    {
        NS_LOG_WARN("Not all UEs have generated CQI feedback.");
    }

    if (!cqiTraces.empty())
    {
        rank = averageRiForAllUes / cqiTraces.size();
        mcs = averageMcsForAllUes / cqiTraces.size();
    }
    else
    {
        rank = 1;
        mcs = 0;
    }

    double flowDuration = (simTime - udpAppStartTime).GetSeconds();
    for (std::map<FlowId, FlowMonitor::FlowStats>::const_iterator i = stats.begin();
         i != stats.end();
         ++i)
    {
        Ipv4FlowClassifier::FiveTuple t = classifier->FindFlow(i->first);
        std::stringstream protoStream;
        protoStream << (uint16_t)t.protocol;
        if (t.protocol == 6)
        {
            protoStream.str("TCP");
        }
        if (t.protocol == 17)
        {
            protoStream.str("UDP");
        }
        outFile << "Flow " << i->first << " (" << t.sourceAddress << ":" << t.sourcePort << " -> "
                << t.destinationAddress << ":" << t.destinationPort << ") proto "
                << protoStream.str() << "\n";
        outFile << "  Tx Packets: " << i->second.txPackets << "\n";
        outFile << "  Tx Bytes:   " << i->second.txBytes << "\n";
        outFile << "  TxOffered:  " << i->second.txBytes * 8.0 / flowDuration / 1000.0 / 1000.0
                << " Mbps\n";
        outFile << "  Rx Bytes:   " << i->second.rxBytes << "\n";
        if (i->second.rxPackets > 0)
        {
            // Measure the duration of the flow from receiver's perspective
            averageFlowThroughput += i->second.rxBytes * 8.0 / flowDuration / 1000 / 1000;
            averageFlowDelay += 1000 * i->second.delaySum.GetSeconds() / i->second.rxPackets;

            double thr = i->second.rxBytes * 8.0 / flowDuration / 1000 / 1000;
            double delay = 1000 * i->second.delaySum.GetSeconds() / i->second.rxPackets;
            double jitter = 1000 * i->second.jitterSum.GetSeconds() / i->second.rxPackets;
            double packetLoss = 1 - ((double)(i->second.rxPackets) / (double)i->second.txPackets);

            outFile << "  Throughput: " << thr << " Mbps\n";
            outFile << "  Mean delay:  " << delay << " ms\n";
            outFile << "  Mean jitter:  " << jitter << " ms\n";
        }
        else
        {
            outFile << "  Throughput:  0 Mbps\n";
            outFile << "  Mean delay:  0 ms\n";
            outFile << "  Mean jitter: 0 ms\n";
        }
        outFile << "  Rx Packets: " << i->second.rxPackets << "\n";
    }

    outFile << "\n\n  Mean flow throughput: " << averageFlowThroughput / stats.size() << "\n";
    outFile << "  Mean flow delay: " << averageFlowDelay / stats.size() << "\n";
    outFile << " Mean rank: " << rank << "\n";
    outFile << " Mean MCS: " << mcs << "\n";

    outFile.close();
    std::ifstream f(filename.c_str());
    if (f.is_open())
    {
        std::cout << f.rdbuf();
    }

    Simulator::Destroy();

    CloseDlSinrFiles();
    return 0;
}
