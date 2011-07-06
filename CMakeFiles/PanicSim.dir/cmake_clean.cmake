FILE(REMOVE_RECURSE
  "./PanicSim_generated_PanicSimMain.cu.o"
  "./PanicSim_generated_hostFunc.cu.o"
  "./PanicSim_generated_kernels.cu.o"
  "./PanicSim_generated_Xlib_mod.cu.o"
  "PanicSim.pdb"
  "PanicSim"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/PanicSim.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
