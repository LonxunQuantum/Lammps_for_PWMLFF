#ifndef LMP_GIT_VERSION_H
#define LMP_GIT_VERSION_H
bool LAMMPS_NS::LAMMPS::has_git_info() { return true; }
const char *LAMMPS_NS::LAMMPS::git_commit() { return "2ef1727a1c5bbeda5778a7f6bf882d7acff493c3"; }
const char *LAMMPS_NS::LAMMPS::git_branch() { return "libtorch"; }
const char *LAMMPS_NS::LAMMPS::git_descriptor() { return ""; }
#endif
