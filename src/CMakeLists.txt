add_library(rsm INTERFACE)

target_compile_features(rsm INTERFACE cxx_std_11)
target_include_directories(rsm INTERFACE .)

target_sources(rsm INTERFACE
  rsm/rsm.hpp
)
