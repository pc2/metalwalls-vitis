# ===============
# Host/test rules
# ===============

CXX = mpic++
CXXFLAGS = -Wall -O3 --std=c++17 -I$(XILINX_XRT)/include -Isrc/ -pthread -fopenmp -I/opt/software/pc2/EB-SW/software/VTK/9.2.2-foss-2022a/include/vtk-9.2
VTK_LIBS = -lvtkCommonCore-9.2 -lvtkCommonDataModel-9.2 -lvtkIOXML-9.2 -lvtksys-9.2
LDFLAGS = -Wall -O3 --std=c++17 -L$(XILINX_XRT)/lib -lOpenCL -pthread -fopenmp $(VTK_LIBS)

COMMON_OBJECTS = build/host/xcl2.o build/host/ComputeUnit.o build/host/Accelerator.o build/host/Experiment.o
COMMON_HEADERS = src/host/xcl2.hpp src/host/Accelerator.hpp src/host/ComputeUnit.hpp src/host/Experiment.hpp

HOST_EXE = build/metalwalls
HOST_OBJECTS = $(COMMON_OBJECTS) build/host/Visualisation.o build/host/main.o
HOST_HEADERS = $(COMMON_HEADERS) src/host/Visualisation.hpp src/host/json.hpp

TEST_EXE = build/unit_tests
TEST_OBJECTS = $(COMMON_OBJECTS) build/test/catch_amalgamated.o build/test/cpu_reference.o build/test/test_k0_acc.o build/test/test_lr_acc.o build/test/test_sr_acc.o build/test/test_data.o build/test/test_distribution_algorithm.o
TEST_HEADERS = $(COMMON_HEADERS) src/test/catch_amalgamated.hpp src/test/KernelRun.hpp

build/host/%.o: src/host/%.cpp $(HOST_HEADERS) Makefile
	mkdir -p build/host/ && $(CXX) $(CXXFLAGS) -c -o $@ $<

build/test/%.o: src/test/%.cpp $(TEST_HEADERS) Makefile
	mkdir -p build/test/ && $(CXX) $(CXXFLAGS) -c -o $@ $<

$(HOST_EXE): $(HOST_OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^

.PHONY: host
host: $(HOST_EXE)

$(TEST_EXE): $(TEST_OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^

.PHONY: test
test: $(TEST_EXE)

# ============
# Device rules
# ============

ifndef $(TARGET)
TARGET = hw
endif

AURORA_KERNELS = AuroraFlow/aurora_flow_0.xo AuroraFlow/aurora_flow_1.xo

.PHONY: aurora
aurora:
	$(MAKE) -C AuroraFlow

VPP = v++
VPP_FLAGS = --target $(TARGET) --platform $(PLATFORM) --save-temps --temp_dir build/.vpptemp
ifeq ($(TARGET),$(filter $(TARGET),sw_emu hw_emu))
	VPP_FLAGS += -DEMULATION
endif

K0_KERNELS= build/k0_acc/kernel_k0_$(TARGET).xo
K0_CONFIG = src/k0_acc/link.cfg

build/k0_acc/%_$(TARGET).xo: src/k0_acc/%.cpp Makefile
	$(VPP) -c $(VPP_FLAGS) -k $* -o $@ '$<'

build/k0_acc_$(TARGET).xclbin: $(K0_KERNELS) $(K0_CONFIG) Makefile
	$(VPP) -l $(VPP_FLAGS) --config $(K0_CONFIG) -o $@ $(K0_KERNELS)

.PHONY: k0_acc
k0_acc: build/k0_acc_$(TARGET).xclbin

LR_KERNELS = build/lr_acc/kernel_lr_$(TARGET).xo
LR_CONFIG = src/lr_acc/link.cfg

build/lr_acc/%_$(TARGET).xo: src/lr_acc/%.cpp Makefile
	$(VPP) -c $(VPP_FLAGS) -k $* -o $@ '$<'

build/lr_acc_$(TARGET).xclbin: $(LR_KERNELS) $(LR_CONFIG) Makefile
	$(VPP) -l $(VPP_FLAGS) --config $(LR_CONFIG) -o $@ $(LR_KERNELS)

.PHONY: lr_acc
lr_acc: build/lr_acc_$(TARGET).xclbin


SR_KERNELS = build/sr_acc/kernel_sr_$(TARGET).xo
SR_CONFIG = src/sr_acc/link.cfg

build/sr_acc/%_$(TARGET).xo: src/sr_acc/%.cpp Makefile
	$(VPP) -c $(VPP_FLAGS) -k $* -o $@ '$<'

build/sr_acc_$(TARGET).xclbin: $(SR_KERNELS) $(SR_CONFIG) Makefile
	$(VPP) -l $(VPP_FLAGS) --config $(SR_CONFIG) -o $@ $(SR_KERNELS)

.PHONY: sr_acc
sr_acc: build/sr_acc_$(TARGET).xclbin

ALL_KERNELS = $(K0_KERNELS) $(LR_KERNELS) $(SR_KERNELS) 
ALL_ACCELERATORS = build/k0_acc_$(TARGET).xclbin build/lr_acc_$(TARGET).xclbin build/sr_acc_$(TARGET).xclbin

.PHONY: reports
reports: $(ALL_KERNELS)

.PHONY: accelerators
accelerators: $(ALL_ACCELERATORS)

# ===============
# Execution rules
# ===============

ifndef $(INPUT)
INPUT = ./experiments/supercap/input_f90.dat
endif

.PHONY: run_app
run_app: $(HOST_EXE)
ifeq ($(TARGET),$(filter $(TARGET),sw_emu hw_emu))
	MW_BINARY_DIR=./build/ XCL_EMULATION_MODE=$(TARGET) mpiexec -n 4 ./$(HOST_EXE) $(INPUT)
else
	MW_BINARY_DIR=./build/ mpiexec -n 4 ./$(HOST_EXE) $(INPUT)
endif

.PHONY: run_test
run_test: $(TEST_EXE)
ifeq ($(TARGET),$(filter $(TARGET),sw_emu hw_emu))
	MW_BINARY_DIR=./build/ XCL_EMULATION_MODE=$(TARGET) ./$(TEST_EXE) $(TEST_PARAMS)
else
	MW_BINARY_DIR=./build/ ./$(TEST_EXE) $(TEST_PARAMS)
endif

.PHONY: clean
clean:
	git clean -Xdf
	
