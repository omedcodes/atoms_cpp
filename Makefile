ifeq ($(OS),Windows_NT)
	SHELL := cmd.exe
	CXX_TEST := $(shell clang++ --version 2>nul)
	ifneq ($(CXX_TEST),)
		CXX = clang++
	else
		CXX = g++
	endif
	EXT = .exe
	BIN = bin\\
	CXXFLAGS = -O2 -std=c++17 -Wall -DGLEW_STATIC
	LDFLAGS = -static -static-libgcc -static-libstdc++ -s -lglfw3 -lglew32 -lopengl32 -lgdi32 -luser32 -lshell32
	MKDIR = if not exist bin mkdir bin
	RM = del /Q /F
else
	CXX_TEST := $(shell clang++ --version 2>/dev/null)
	ifneq ($(CXX_TEST),)
		CXX = clang++
	else
		CXX = g++
	endif
	EXT =
	BIN = bin/
	CXXFLAGS = -O2 -flto -std=c++17 -Wall -DGLEW_STATIC
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		LDFLAGS = -lglfw -lGLEW -lGL -lm
	endif
	ifeq ($(UNAME_S),Darwin)
		LDFLAGS = -lglfw -lGLEW -framework OpenGL
	endif
	MKDIR = mkdir -p bin
	RM = rm -rf bin/
endif

TARGET1 = $(BIN)atom_raytracer$(EXT)
TARGET2 = $(BIN)atome_wave2d$(EXT)
SRC1 = src/atom_raytracer.cpp
SRC2 = src/atom_wave2d.cpp

all: $(TARGET1) $(TARGET2)

$(TARGET1): $(SRC1)
	@$(MKDIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(TARGET2): $(SRC2)
	@$(MKDIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	$(RM) $(BIN)*$(EXT)