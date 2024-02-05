#ifndef LMP_GIT_VERSION_H
#define LMP_GIT_VERSION_H
bool LAMMPS_NS::LAMMPS::has_git_info() { return true; }
const char *LAMMPS_NS::LAMMPS::git_commit() { return "5a198eee55984c8fcef36b5a242f357e779b3214"; }
const char *LAMMPS_NS::LAMMPS::git_branch() { return "libtorch"; }
const char *LAMMPS_NS::LAMMPS::git_descriptor() { return ""; }
#endif
