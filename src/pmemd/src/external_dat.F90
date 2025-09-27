! Another module needed to avoid cyclic dependencies

module external_dat_mod

  implicit none

  integer, parameter            :: max_fn_len = 256

  public

  character(max_fn_len), save :: extprog, json

end module external_dat_mod
