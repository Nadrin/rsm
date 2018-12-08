add_library(rsm INTERFACE)

target_compile_features(rsm INTERFACE cxx_std_14)
target_include_directories(rsm INTERFACE .)

target_sources(rsm INTERFACE
  rsm/rsm.hpp
  rsm/init.hpp
  rsm/lds.hpp
  rsm/next.hpp
  rsm/options.hpp
  rsm/range.hpp
  rsm/utils.hpp
  rsm/detail/common.hpp
  rsm/detail/ldsperm.hpp
  rsm/detail/memory.hpp
  rsm/detail/primes.hpp
  rsm/detail/rijumptable.inl
  rsm/distributions/disk.hpp
  rsm/distributions/hemisphere.hpp
  rsm/distributions/ncube.hpp
  rsm/distributions/sphere.hpp
  rsm/generators/pcg32.hpp
  rsm/generators/splitmix64.hpp
  rsm/generators/stlcompat.hpp
  rsm/generators/xoroshiro128p.hpp
  rsm/generators/xoroshiro64s.hpp
  rsm/samplers/halton.hpp
  rsm/samplers/hammersley.hpp
  rsm/samplers/lhs.hpp
  rsm/samplers/random.hpp
  rsm/samplers/stratified.hpp
  rsm/montecarlo/heuristics.hpp
)
