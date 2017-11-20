# Scaffold_graph
Creates scaffold graph (bioinformatic thing)

How to use it:
1. Create static-library 'libyaml-cpp.a' of YAML-0.3 (https://github.com/jbeder/yaml-cpp/releases/tag/release-0.3.0) and copy it in the root of the project.
2. Copy folder 'yaml-cpp/include/yaml-cpp' to the root.
3. Create the YAML-file (example.yaml you can find in the repository) 
4. CMake it!
5. ???
6. Profit! 

Output is graph in .dot format and optionally (depends on your YAML) some helpful for next runs files of mine-defined format (.his and .mtrx).
Also apart from YAML-file you have to have at least one pair (left and right) of .sam files which I call "library" (it is clear from example of YAML). They are huge, so I don't include them in my repository for an example.
