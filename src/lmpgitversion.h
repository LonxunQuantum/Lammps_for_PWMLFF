#ifndef LMP_GIT_VERSION_H
#define LMP_GIT_VERSION_H
bool LAMMPS_NS::LAMMPS::has_git_info() { return true; }
const char *LAMMPS_NS::LAMMPS::git_commit() { return "4a1ef507be0b6505d4795e6e7aef4ec1e7458e50"; }
const char *LAMMPS_NS::LAMMPS::git_branch() { return "libtorch"; }
const char *LAMMPS_NS::LAMMPS::git_descriptor() { return ""; }
#endif
