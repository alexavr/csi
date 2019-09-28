# Description
Currently computes eulerian criteria only: $Q$,$\Delta$ and $\lambda2$.
Working on lagrangian approach.

# Tech notes

##### How to add new criterion

1. Add new logical variable into `namelist.csi`

   Try to keep the convention of variable name: *{criterion_name}_criterion*, by default all criteria are set to be *true*. Keep it this way.

2. Describe it in `module_globals.f90`

3. Add calculation subroutine into `module_csi.f90`

   The stress tensor calculated before so you don't need to do it again. Be aware all computation performed in real(8), but the result truncated into real(4) in the sake of storage. The subroutine has to have output (see another subroutine as an example).

4. Recompile, run and enjoy

