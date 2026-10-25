# MU-MIMO with Type-II Codebook Precoding and User Grouping in ns-3 5G-LENA

<!-- TODO: replace with the final paper title, venue, and badges -->
[![ns-3](https://img.shields.io/badge/ns--3-3.xx-blue)](https://www.nsnam.org/)
[![5G-LENA](https://img.shields.io/badge/5G--LENA-v4.x-green)](https://cttc-lena.gitlab.io/nr/html/)
[![License](https://img.shields.io/badge/license-GPL--2.0-lightgrey)](LICENSE)

This repository contains the simulation code accompanying the paper:

> **\<MU-MIMO with Type-II Codebook Precoding and User Grouping in ns-3 5G-LENA\>**
> Amirhassan Babazadeh Darabi, \<Co-authors\>, Marina Petrova, Ljiljana Simić
> *\<MobiWac, 2026\>*
> Chair of Mobile Communications and Computing (MCC), RWTH Aachen University

It extends [5G-LENA](https://gitlab.com/cttc-lena/nr) with multi-user MIMO (MU-MIMO) spatial-division multiple access (SDMA) scheduling based on **CPOG — Codebook-based Precoder Orthogonality Grouping**. CPOG groups users for co-scheduling using only the precoding matrix indicators reported through Type-II codebook CSI feedback, and serves each group with zero-forcing (ZF) precoding, avoiding the need for full channel knowledge at the gNB.

---

## Table of Contents

- [Features](#features)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Running the Simulation](#running-the-simulation)
- [Parameters](#parameters)
- [Output](#output)
- [Reproducing the Paper Results](#reproducing-the-paper-results)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Features

- **Type-II codebook** (`NrCbTypeTwo`) and corresponding **PMI search** (`NrPmSearchTypeTwo`) for high-resolution CSI feedback
- **TDMA + SDMA MAC scheduler** (`NrMacSchedulerTdmaSdma`, round-robin variant `NrMacSchedulerTdmaSdmaRR`) supporting multiple co-scheduled users per slot
- **User grouping methods**
  - CPOG (Type-I and Type-II codebook-based precoder orthogonality)
  - Full-CSI SVD baseline
  - Volume-ratio baseline
  - K-means clustering
- **ZF precoding** computed from Type-II codebook feedback
- **Post-precoding SINR computation** for co-scheduled users
- Per-RNTI dynamic BWP routing via a custom `BwpManagerAlgorithm` subclass
- Python analysis scripts for SINR CDFs, co-scheduling degree, per-UE throughput, and scalability plots

---

## Repository Structure

<!-- TODO: adjust to the actual layout of the repository -->
```
.
├── examples/
│   └── mcc-nr-mu-mimo.cc          # Main simulation script
├── model/                         # Modified / new 5G-LENA model files
│   ├── nr-cb-type-two.{h,cc}
│   ├── nr-pm-search-type-two.{h,cc}
│   ├── nr-mac-scheduler-tdma-sdma.{h,cc}
│   └── nr-mac-scheduler-tdma-sdma-rr.{h,cc}
├── scripts/
│   ├── analysis/                  # RxPacketTrace post-processing and plotting
│   └── slurm/                     # Batch scripts for parameter sweeps
└── README.md
```

---

## Requirements

| Component | Version |
|-----------|---------|
| ns-3      | `ns-3.xx` |
| 5G-LENA   | `5g-lena-v4.x.y` |
| Compiler  | GCC ≥ 11 or Clang ≥ 14 |
| CMake     | ≥ 3.13 |
| Eigen3    | Required (MIMO / codebook computations) |
| Python    | ≥ 3.8 (analysis scripts only) |

Python packages for analysis: `numpy`, `pandas`, `matplotlib`, `scipy`.

---

## Installation

```bash
# 1. Clone ns-3
git clone https://gitlab.com/nsnam/ns-3-dev.git
cd ns-3-dev
git checkout ns-3.xx

# 2. Clone this repository as the nr module
cd contrib
git clone https://github.com/<user>/<repo>.git nr
cd ..

# 3. Install Eigen (Ubuntu/Debian)
sudo apt install libeigen3-dev

# 4. Configure and build
./ns3 configure --enable-examples --enable-tests -d optimized
./ns3 build
```

> **Tip:** Use the `optimized` profile for parameter sweeps; switch to `-d debug` when developing.

---

## Running the Simulation

The reference configuration used in the paper:

```bash
./ns3 run "mcc-nr-mu-mimo.cc \
  --pmSearchMethod=ns3::NrPmSearchTypeTwo \
  --fullSearchCb=ns3::NrCbTypeTwo \
  --rankLimit=1 \
  --subbandSize=8 \
  --numRowsGnb=8 --numColumnsGnb=4 \
  --numVPortsGnb=4 --numHPortsGnb=4 --xPolGnb=true \
  --numRowsUe=4 --numColumnsUe=2 \
  --numVPortsUe=2 --numHPortsUe=2 --xPolUe=true \
  --packetInterval=160us \
  --simTime=1 \
  --bandwidth=10e6 \
  --scheduler=ns3::NrMacSchedulerTdmaSdmaRR \
  --beamformingMethod=ns3::DirectPathBeamforming \
  --Users=4 --maxSdmaUsers=4 --MaxStreams=4 \
  --runNumber=42 \
  --threshold=3e-1 \
  --csiFlags=6 \
  --ns3::NrMacSchedulerTdmaSdma::FullCsiBaseline=false"
```

To run the **Full-CSI SVD baseline** instead of CPOG, set:

```bash
--ns3::NrMacSchedulerTdmaSdma::FullCsiBaseline=true
```

<!-- TODO: add the attribute/flag used to select the volume-ratio and k-means baselines -->

---

## Parameters

### CSI and Codebook

| Parameter | Value | Description |
|-----------|-------|-------------|
| `pmSearchMethod` | `ns3::NrPmSearchTypeTwo` | PMI search algorithm used by the UE |
| `fullSearchCb` | `ns3::NrCbTypeTwo` | Codebook used for the PMI search (Type-II) |
| `rankLimit` | `1` | Maximum transmission rank per UE |
| `subbandSize` | `8` | Subband size in PRBs for subband PMI/CQI reporting |
| `csiFlags` | `6` | CSI feedback flags bitmask (`CQI_PDSCH_MIMO` \| `CQI_CSI_RS`) |

### Antenna Configuration

| Parameter | gNB | UE |
|-----------|-----|----|
| Rows × Columns | 8 × 4 | 4 × 2 |
| Vertical ports | 4 | 2 |
| Horizontal ports | 4 | 2 |
| Dual-polarized | Yes | Yes |
| **Total antenna ports** | **32** | **8** |

The number of ports is `numVPorts × numHPorts × (xPol ? 2 : 1)`.

### Scheduling and MU-MIMO

| Parameter | Value | Description |
|-----------|-------|-------------|
| `scheduler` | `ns3::NrMacSchedulerTdmaSdmaRR` | TDMA+SDMA round-robin scheduler |
| `beamformingMethod` | `ns3::DirectPathBeamforming` | Initial analog beamforming method |
| `Users` | `4` | Number of UEs in the scenario |
| `maxSdmaUsers` | `4` | Maximum number of users co-scheduled in one SDMA group |
| `MaxStreams` | `4` | Maximum total spatial streams per slot |
| `threshold` | `0.03` | Orthogonality threshold for CPOG user grouping |
| `NrMacSchedulerTdmaSdma::FullCsiBaseline` | `false` | `false` = CPOG (codebook-based), `true` = Full-CSI SVD baseline |

### Traffic and Simulation

| Parameter | Value | Description |
|-----------|-------|-------------|
| `packetInterval` | `160us` | UDP packet inter-arrival time (near full-buffer load) |
| `bandwidth` | `10e6` | Channel bandwidth in Hz (10 MHz) |
| `simTime` | `1` | Simulation duration in seconds |
| `runNumber` | `42` | RNG run number for reproducibility |

---

## Output

The simulation writes the standard 5G-LENA trace files to the working directory, most importantly:

- `RxPacketTrace.csv` — per-TB reception records (RNTI, SINR, MCS, TB size, corruption status), used for all SINR and throughput analysis
- <!-- TODO: list any additional custom traces, e.g. grouping decisions or co-scheduling logs -->

---

## Reproducing the Paper Results

### 1. Parameter sweep

The scalability results sweep the number of users from K = 1 to K = 16 over multiple independent runs. Example local loop:

```bash
for K in 1 2 4 8 12 16; do
  for RUN in $(seq 1 50); do
    ./ns3 run "mcc-nr-mu-mimo.cc <reference flags> --Users=$K --runNumber=$RUN" \
      --cwd="results/K${K}_run${RUN}"
  done
done
```


### 2. Analysis

```bash
cd scripts/analysis
python3 plot_sinr_cdf.py        --input ../../results   # SINR CDFs
python3 plot_coscheduling.py    --input ../../results   # Co-scheduling degree
python3 plot_throughput_box.py  --input ../../results   # Per-UE throughput box plots
python3 plot_scalability.py     --input ../../results   # Scalability vs. K
```

<!-- TODO: match script names and arguments to the actual files -->

---

## Citation

If you use this code in your research, please cite:


---

## License

This project builds on ns-3 and 5G-LENA and is distributed under the **GNU General Public License v2.0**, consistent with the upstream projects. See [LICENSE](LICENSE) for details.

---

## Contact

**Amirhassan Babazadeh Darabi**
Chair of Mobile Communications and Computing (MCC), RWTH Aachen University


Issues and pull requests are welcome.
