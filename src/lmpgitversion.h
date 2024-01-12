#ifndef LMP_GIT_VERSION_H
#define LMP_GIT_VERSION_H
bool LAMMPS_NS::LAMMPS::has_git_info() { return true; }
const char *LAMMPS_NS::LAMMPS::git_commit() { return "6b7fd2f59d5359245e51b600d8c91ce48f9a85df"; }
const char *LAMMPS_NS::LAMMPS::git_branch() { return "libtorch"; }
const char *LAMMPS_NS::LAMMPS::git_descriptor() { return ""; }
#endif
