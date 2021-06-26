#
# Copyright 2019 Xilinx, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
############################## Help Section ##############################
.PHONY: help

help::
	$(ECHO) "Makefile Usage:"
	$(ECHO) "  make all TARGET=<sw_emu/hw_emu/hw> DEVICE=<FPGA platform> HOST_ARCH=<aarch32/aarch64/x86>"
	$(ECHO) "      Command to generate the design for specified Target and Shell."
	$(ECHO) "      By default, HOST_ARCH=x86. HOST_ARCH is required for SoC shells"
	$(ECHO) ""
	$(ECHO) "  make clean "
	$(ECHO) "      Command to remove the generated non-hardware files."
	$(ECHO) ""
	$(ECHO) "  make cleanall"
	$(ECHO) "      Command to remove all the generated files."
	$(ECHO) ""
	$(ECHO) "  make sd_card TARGET=<sw_emu/hw_emu/hw> DEVICE=<FPGA platform> HOST_ARCH=<aarch32/aarch64/x86>"
	$(ECHO) "      Command to prepare sd_card files."
	$(ECHO) "      By default, HOST_ARCH=x86. HOST_ARCH is required for SoC shells"
	$(ECHO) ""
	$(ECHO) "  make run TARGET=<sw_emu/hw_emu/hw> DEVICE=<FPGA platform> HOST_ARCH=<aarch32/aarch64/x86>"
	$(ECHO) "      Command to run application in emulation."
	$(ECHO) "      By default, HOST_ARCH=x86. HOST_ARCH required for SoC shells"
	$(ECHO) ""
	$(ECHO) "  make xclbin TARGET=<sw_emu/hw_emu/hw> DEVICE=<FPGA platform> HOST_ARCH=<aarch32/aarch64/x86>"
	$(ECHO) "      Command to build xclbin application."
	$(ECHO) "      By default, HOST_ARCH=x86. HOST_ARCH is required for SoC shells"
	$(ECHO) ""
	$(ECHO) "  make host DEVICE=<FPGA platform> HOST_ARCH=<aarch32/aarch64/x86>"
	$(ECHO) "      Command to build host application."
	$(ECHO) "      By default, HOST_ARCH=x86. HOST_ARCH is required for SoC shells"
	$(ECHO) ""
	$(ECHO) "  NOTE: For SoC shells, ENV variable SYSROOT needs to be set."
	$(ECHO) ""

############################## Setting up Project Variables ##############################
MK_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
XF_PROJ_ROOT ?= $(shell bash -c 'export MK_PATH=$(MK_PATH); echo $${MK_PATH%L2/tests/gelinearsolver/*}')
CUR_DIR := $(patsubst %/,%,$(dir $(MK_PATH)))
XFLIB_DIR = $(XF_PROJ_ROOT)

TARGET ?= hw
HOST_ARCH := x86
SYSROOT := ${SYSROOT}
DEVICE ?= xilinx_u50_gen3x16_xdma_201920_3

RUNS ?= 1
SEED ?= 12

NUM_CUS ?= 36
MAX_SIZE ?= 500

include ./utils.mk

XDEVICE := $(call device2xsa, $(DEVICE))
TEMP_DIR := _x_temp.$(TARGET).$(XDEVICE)
TEMP_REPORT_DIR := $(CUR_DIR)/reports/_x.$(TARGET).$(XDEVICE)
BUILD_DIR := build_dir.$(TARGET).$(XDEVICE)
BUILD_REPORT_DIR := $(CUR_DIR)/reports/_build.$(TARGET).$(XDEVICE)
EMCONFIG_DIR := $(BUILD_DIR)

# Setting tools
VPP := v++
SDCARD := sd_card
EMU_DIR := $(SDCARD)/data/emulation

############################## Setting up Host Variables ##############################
#Include Required Host Source Files
#HOST_SRCS += test_gelinearsolver_interface.cpp
#HOST_SRCS += test_gelinearsolver.cpp
HOST_SRCS += TutorialCpp_main.cpp
#HOST_SRCS += $(XFLIB_DIR)/ext/xcl2/xcl2.cpp
#HOST_SRCS += IpVitisSolverInterface.cpp
#HOST_SRCS += MyNLP.cpp
#HOST_SRCS += hs071_nlp.cpp
#HOST_SRCS += IpAlgBuilder.cpp
HOST_SRCS += TutorialCpp_nlp.cpp

#CXXFLAGS += -I$(XFLIB_DIR)/L2/tests/gelinearsolver
#CXXFLAGS += -I$(XFLIB_DIR)/ext/xcl2
#CXXFLAGS += -I$(XFLIB_DIR)/ext/MatrixGen
#CXXFLAGS += -I$(XFLIB_DIR)/L2/tests/gelinearsolver/Ipopt-releases-3.13.3/src/Algorithm/LinearSolvers
#CXXFLAGS += -I$(XFLIB_DIR)/L2/tests/gelinearsolver/Ipopt-releases-3.13.3/src/Algorithm
#CXXFLAGS += -I$(XFLIB_DIR)/L2/tests/gelinearsolver/Ipopt-releases-3.13.3/src/LinAlg
#CXXFLAGS += -I$(XFLIB_DIR)/L2/tests/gelinearsolver/Ipopt-releases-3.13.3/src/contrib/CGPenalty
CXXFLAGS += -I/usr/local/include/coin-or

#INCL = `PKG_CONFIG_PATH=/usr/local/lib/pkgconfig: pkg-config --cflags ipopt`
#LIBS = `PKG_CONFIG_PATH=/usr/local/lib/pkgconfig: pkg-config --libs ipopt`
#CXXLINKFLAGS =  -Wl,--rpath -Wl,/usr/local/lib



# Host compiler global settings
CXXFLAGS += -I$(XILINX_XRT)/include -I$(XILINX_VIVADO)/include -std=c++14  -Wall -Wno-unknown-pragmas -Wno-unused-label -pthread 
CXXLINKFLAGS += -L$(XILINX_XRT)/lib -lxilinxopencl -lpthread -lrt -Wno-unused-label -Wno-narrowing -DVERBOSE 
CXXLINKFLAGS += -L$(XILINX_VIVADO)/lnx64/tools/fpo_v7_0 -Wl,--as-needed -lgmp -lmpfr -lIp_floating_point_v7_0_bitacc_cmodel
CXXLINKFLAGS += -L/usr/local/lib -lipopt -lvitislibrary
CXXFLAGS += -fmessage-length=0 -O2 -DNDEBUG -g
#CXXFLAGS += -g -O0

#EXE_FILE := test_gelinearsolver_interface.exe
#EXE_FILE := test_gelinearsolver.exe
EXE_FILE := TutorialCpp_main.exe
SOC_HOST_ARGS :=  -xclbin $(BUILD_DIR)/kernel_gelinearsolver.xclbin
HOST_ARGS += -xclbin $(BUILD_DIR)/kernel_gelinearsolver.xclbin
HOST_ARGS += -M $(ROW_SIZE)
HOST_ARGS += -N $(COLUMN_SIZE)

ifneq ($(HOST_ARCH), x86)
	LDFLAGS += --sysroot=$(SYSROOT)
endif
############################## Setting up Kernel Variables ##############################
# Kernel compiler global settings
VPP_FLAGS += -t $(TARGET) --platform $(XPLATFORM) --save-temps
LDCLFLAGS += --jobs 8
ifneq (,$(shell echo $(XPLATFORM) | awk '/u280/'))
# U280
VPP_FLAGS += --config conn_u280.ini
CXXFLAGS += -DUSE_HBM
else ifneq (,$(shell echo $(XPLATFORM) | awk '/u250/'))
# U250
VPP_FLAGS += --config conn_u250.ini
CXXFLAGS += -DUSE_DDR
else ifneq (,$(shell echo $(XPLATFORM) | awk '/u200/'))
# U200
VPP_FLAGS += --config conn_u200.ini
CXXFLAGS += -DUSE_DDR
else ifneq (,$(shell echo $(XPLATFORM) | awk '/u50/'))
# U50
VPP_FLAGS += --config conn_u50.ini
CXXFLAGS += -DUSE_HBM
else ifneq (,$(XPLATFORM))
$(warning Unsupported platform $(XPLATFORM))
endif

ifeq ($(TARGET),sw_emu)
VPP_CFLAGS += -D VIVADO_HLS_SIM
endif
ifeq ($(TARGET),hw_emu)
VPP_CFLAGS += -D VIVADO_HLS_SIM
endif
VPP_FLAGS += -I$(XFLIB_DIR)/L2/tests/gelinearsolver
VPP_FLAGS += -I$(XFLIB_DIR)/L2/include
VPP_FLAGS += -I$(XFLIB_DIR)/ext

VPP_FLAGS += --optimize 3
VPP_FLAGS += --kernel_frequency 500
VPP_FLAGS += -D NCU_IN=$(NUM_CUS)
VPP_FLAGS += -D MAXN_IN=$(MAX_SIZE)


############################## Declaring Binary Containers ##############################
BINARY_CONTAINERS += $(BUILD_DIR)/kernel_gelinearsolver.xclbin
BINARY_CONTAINER_kernel_gelinearsolver_OBJS += $(TEMP_DIR)/kernel_gelinearsolver_0.xo

############################## Setting Targets ##############################
CP = cp -rf

.PHONY: all clean cleanall docs emconfig
all: check_vpp check_platform | $(EXE_FILE) $(BINARY_CONTAINERS) emconfig sd_card

.PHONY: host
host: $(EXE_FILE) | check_xrt

.PHONY: xclbin
xclbin: check_vpp | $(BINARY_CONTAINERS)

############################## Setting Rules for Binary Containers (Building Kernels) ##############################
$(TEMP_DIR)/kernel_gelinearsolver_0.xo: $(XFLIB_DIR)/L2/tests/gelinearsolver/kernel_gelinearsolver.cpp
	$(ECHO) "Compiling Kernel: kernel_gelinearsolver_0"
	mkdir -p $(TEMP_DIR)
	$(VPP) $(kernel_gelinearsolver_0_VPP_FLAGS) $(VPP_FLAGS) --temp_dir $(TEMP_DIR) --report_dir $(TEMP_REPORT_DIR) -c -k kernel_gelinearsolver_0 -I'$(<D)' -o $@ $^

$(BUILD_DIR)/kernel_gelinearsolver.xclbin: $(BINARY_CONTAINER_kernel_gelinearsolver_OBJS)
	mkdir -p $(BUILD_DIR)
	$(VPP) $(VPP_FLAGS) --temp_dir $(BUILD_DIR) --report_dir $(BUILD_REPORT_DIR)/kernel_gelinearsolver -l $(LDCLFLAGS) $(LDCLFLAGS_kernel_gelinearsolver) -o'$@' $(+)

############################## Setting Rules for Host (Building Host Executable) ##############################
prefix=/usr/local
exec_prefix=${prefix}
 
$(EXE_FILE): $(HOST_SRCS)
	$(CXX) -o $@ $^ $(CXXLINKFLAGS) $(INCL) $(CXXFLAGS)

emconfig:$(EMCONFIG_DIR)/emconfig.json
$(EMCONFIG_DIR)/emconfig.json:
	emconfigutil --platform $(XPLATFORM) --od $(EMCONFIG_DIR)

############################## Setting Essential Checks and Running Rules ##############################
run: all
ifeq ($(findstring zc*, $(DEVICE)), zc*)
$(error This example is not supported for $(DEVICE))
endif

ifeq ($(TARGET),$(filter $(TARGET),sw_emu hw_emu))
ifeq ($(HOST_ARCH), x86)
	$(CP) $(EMCONFIG_DIR)/emconfig.json .
	XCL_EMULATION_MODE=$(TARGET) ./$(EXE_FILE) $(HOST_ARGS)
else
	mkdir -p $(EMU_DIR)
	$(CP) $(XILINX_VITIS)/data/emulation/unified $(EMU_DIR)
	mkfatimg $(SDCARD) $(SDCARD).img 500000
	launch_emulator -no-reboot -runtime ocl -t $(TARGET) -sd-card-image $(SDCARD).img -device-family $(DEV_FAM)
endif
else
ifeq ($(HOST_ARCH), x86)
	./$(EXE_FILE) $(HOST_ARGS)
else
	$(ECHO) "Please copy the content of sd_card folder and data to an SD Card and run on the board"
endif
endif

############################## Preparing sdcard folder ##############################
sd_card: $(EXE_FILE) $(BINARY_CONTAINERS) emconfig
ifneq ($(HOST_ARCH), x86)
	mkdir -p $(SDCARD)/$(BUILD_DIR)
	mkdir -p $(SDCARD)/data
	$(CP) $(B_NAME)/sw/$(XDEVICE)/boot/generic.readme $(B_NAME)/sw/$(XDEVICE)/xrt/image/* xrt.ini $(EXE_FILE) $(SDCARD)
	$(CP) $(BUILD_DIR)/*.xclbin $(SDCARD)/$(BUILD_DIR)/
ifeq ($(TARGET),$(filter $(TARGET),sw_emu hw_emu))
	$(ECHO) 'cd /mnt/' >> $(SDCARD)/init.sh
	$(ECHO) 'export XILINX_VITIS=$$PWD' >> $(SDCARD)/init.sh
	$(ECHO) 'export XCL_EMULATION_MODE=$(TARGET)' >> $(SDCARD)/init.sh
	$(ECHO) './$(EXE_FILE) $(SOC_HOST_ARGS)' >> $(SDCARD)/init.sh
	$(ECHO) 'reboot' >> $(SDCARD)/init.sh
else
	[ -f $(SDCARD)/BOOT.BIN ] && echo "INFO: BOOT.BIN already exists" || $(CP) $(BUILD_DIR)/sd_card/BOOT.BIN $(SDCARD)/
	$(ECHO) './$(EXE_FILE) $(SOC_HOST_ARGS)' >> $(SDCARD)/init.sh
endif
endif

############################## Cleaning Rules ##############################
cleanh:
	-$(RMDIR) $(EXE_FILE) vitis_* TempConfig system_estimate.xtxt *.rpt .run/
	-$(RMDIR) src/*.ll _xocc_* .Xil dltmp* xmltmp* *.log *.jou *.wcfg *.wdb

cleank:
	-$(RMDIR) $(BUILD_DIR)/*.xclbin _vimage *xclbin.run_summary qemu-memory-_* emulation/ _vimage/ pl* start_simulation.sh *.xclbin
	-$(RMDIR) _x_temp.*/_x.* _x_temp.*/.Xil _x_temp.*/profile_summary.* sample_* 
	-$(RMDIR) _x_temp.*/dltmp* _x_temp.*/kernel_info.dat _x_temp.*/*.log 
	-$(RMDIR) _x_temp.* 

cleanall: cleanh cleank
	-$(RMDIR) $(BUILD_DIR) sd_card* build_dir.* emconfig.json *.html $(TEMP_DIR) *.txt.* $(CUR_DIR)/reports *.csv *.run_summary $(CUR_DIR)/*.raw
	-$(RMDIR) $(XFLIB_DIR)/common/data/*.xe2xd* $(XFLIB_DIR)/common/data/*.orig*


clean: cleanh
