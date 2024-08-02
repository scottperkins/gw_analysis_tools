file(REMOVE_RECURSE
  "../lib/libgwat.a"
  "../lib/libgwat.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/gwat_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
