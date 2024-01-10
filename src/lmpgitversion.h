#ifndef LMP_GIT_VERSION_H
#define LMP_GIT_VERSION_H
bool LAMMPS_NS::LAMMPS::has_git_info() { return true; }
const char *LAMMPS_NS::LAMMPS::git_commit() { return "cf65c0fecf6cc28f366edd2bb3ca8e60e352aeb4"; }
const char *LAMMPS_NS::LAMMPS::git_branch() { return "libtorch"; }
const char *LAMMPS_NS::LAMMPS::git_descriptor() { return ""; }
#endif
