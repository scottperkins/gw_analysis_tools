file(REMOVE_RECURSE
  "../lib/libgwat.pdb"
  "../lib/libgwat.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/gwat.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
