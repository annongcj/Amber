
!*******************************************************************************
!
! Module:  findmask
!
! Description:
!              
! This subroutine takes "atomic expression" loosely using Chimera/Midas syntax
! and decomposes it into series of elementary selections that are connected
! by binary (or unary) operators. The resulting selection is stored in a 
! integer mask array.
!
! Parentheses and logical operators (precedence: ! > & > |) are allowed.
! Expression parsing is done through several intermediate stages: first, 
! the atomic expression is 'tokenized', i.e. 'elementary selections' are 
! enclosed into brackets [..], and basic error checking (e.g. for unknown
! symbols) is done. Second, tokenized expression is converted into postfix 
! notation (or Reverse Polish notation) which gets rid of parentheses
! and defines the order of operations based on operator precedence.
! Finally, postfix notation needs to be evaluated. The resulting selection
! is stored in an integer mask[] array ('0' or '1' for each atom).
! Steps (2) and (3) are done through stack data structure as originally 
! described in one of the Knuth's papers (as well as many other textbooks, 
! e.g. L.R.Nyhoff, S.C.Leestma: Fortran77 for Engineers and Scientists).
!
! The syntax for elementary selections is the following
! :{residue numlist}          e.g. [:1-10] [:1,3,5] [:1-3,5,7-9]
! :{residue namelist}         e.g. [:LYS] [:ARG,ALA,GLY]
! @{atom numlist}             e.g. [@12,17] [@54-85] [@12,54-85,90]
! @{atom namelist}            e.g. [@CA] [@CA,C,O,N,H]
! @%{atom amber typelist}     e.g. [@%CT] [@%N*,N3]
! :* means all residues and @* means all atoms
!                             e.g. [:10@*], or equally [:10]
! 
! A wildcard '=' following a letter (or letters) in atom or residue names
! matches everything that starts with the letter (or letters), e.g.:
! [@C=] matches all atoms that start with C (i.e. effectively all carbons)
! [@H=] matches all types of hydrogens (i.e. H, HH32, HA, etc.)
! Therefore, it can be used to match certain atomic element types, such as
! all hydrogens, or sulfurs, or carbons, etc.
! [:A=] matches all ALA, ASP, ASN and ARG
!
! compound expressions of the following type are also allowed:
! :{residue numlist|namelist}@{atom numlist|namelist|typelist}
! and are processed as (i.e. replaced by two AND'ed expressions):
! :{residue numlist|namelist} & @{atom numlist|namelist|typelist}
! e.g.  :1-10@CA    is equivalent to   :1-10 & @CA
!       :LYS@H=     is equivalent to   :LYS & @H=
!
! more examples:
! :ALA,TRP     ... all alanine and tryptophane residues
! :5,10@CA     ... CA carbon in residues 5 and 10 
! :* & !@H=    ... all non-hydrogen atoms (equivalent to "!@H=")
! @CA,C,O,N,H  ... all backbone atoms
! !@CA,C,O,N,H ... all non-backbone atoms (=sidechains for proteins)
! :1-500@O & !(:WAT | :LYS,ARG)
!              ... all backbone oxygens in residues 1-500 but not in 
!                  water, lysine or arginine residues
!
! Distance operators (<, >) introduce a new kind of 'distance operand',
! which also encodes whether distance cutoff is applied by atoms or by
! residues. The 'distance operand' starts with ':' or '@' followed by
! an integer or floating point number, e.g. [:3.4], [@5.2], [@3].
! Distance operators (<, >) have highest priority, i.e. higher than '!' 
! but still lower than parentheses.
! Here are a few examples of how to use it:
! :1-5<@2.5      ... all *atoms* within 2.5A from residues 1-5
! :1-5 & :ALA<:8 ... intersection of residues 1-5 and all *residues* 
!                    within 8A of residue 7
! (:1-5&:ALA)<:8 ... all *residues* within 8A of all alanines in 
!                    residues 1-5
! Implementation note: Because 'distance operands' are different from
! other operands which always represent elementary masks, it cannot be 
! represented by mask(1..natom) integer array on the stack and has to
! be delt with in a different way (see the code in routine eval() ).
!
! Name selections are case insensitive. Most names in amber topology
! file are in capital letters anyway. The only disadvantage is that
! we cannot distinguish atom names for CA (alpha carbon) and Ca
! (calcium) and some similar rare cases.
!
! Assumptions about residue, atom and element/amber atom type names:
! - atom/residue names are 4 chars long (padded with blanks)
! - names consist of: [a-zA-Z] (letters); [0-9] (digits: in atom names,
!   e.g. HH31, CH3, in amber atom types, e.g. N3, H1, O2, in nucleic
!   acid residues, e.g. DC3, DG5, maybe in other custom made residues);
!   ['] (prime: in nucleic acids, e.g. C5', H5'1); [*] (star: in amber
!   atom types, I think just two of them N*, C*); [+-] (plus and minus:
!   some ions (defined in leap, or custom) will create atom and residue
!   names with + or -, e.g. Cl-, Cs+, K+, Li+, Na+, Rb+).
! - selection of atom elements can be achieved by using '=' wildcard:
!   one (C,N,O) or, less frequently, two character (Na,Cl,Mg,etc.) names
!   effectively represent atomic elements when followed by '='
! - some static buffers during processing of 'maskstr' are set up 
!   and this limits the length of selection string to 256 chars
!
! VH; Dec 30, 2003
! 
! TODO: (things which could be useful to add)
! - LES copy selection (syntax maybe: [#1] [#1-3,7-9], etc.)
!   compound selections should be possible as well, e.g:
!   :35-48#1,3-5 ... copies 1,3,4,5 in residues 35-48
! - Defining selections that could be used in later selections
!   that would be convenient for backbone, sidechains, etc. but
!   maybe implement as {selection name}
! - '?' wildcard: a simple form of '?' in the function isNameMatch
!   should not be difficult but I don't see much use for it
!   '*' is probably not a good idea because it occurs in two amber
!   atom type names (see above)
!*******************************************************************************

module findmask_mod


! Make everything here private except atommask
private
public :: atommask

! Buffer length
integer, parameter :: BUFLEN = 256

contains

!*******************************************************************************
!
! Subroutine:  atommask
!
! Description: 
!     
! Parses maskstr and sets atom selection array. Value of 0 indicates the atom
! is not selected, 1 indicates that it is.
!              
!*******************************************************************************

subroutine atommask(natom,nres,prnlev,igraph,isymbl,ipres,  &
                    lbres,crd,maskstr,mask) 
   implicit none
   
! Passed variables

   integer, intent(in)   :: natom          ! number of atoms
   integer, intent(in)   :: nres           ! number of residues
   integer, intent(in)   :: prnlev         ! print level (0-3)
   integer, intent(in)   :: ipres(nres)    ! residue pointer array

   integer, intent(out)  :: mask(natom)    ! mask selection array

   character (len=4), intent(in)     :: igraph(natom)  ! atom names
   character (len=4), intent(in)     :: isymbl(natom)  ! atom types
   character (len=4), intent(in)     :: lbres(nres)    ! residue labels

   character(len=BUFLEN), intent(inout) :: maskstr       ! mask string

   double precision, intent(in)      :: crd(3*natom)     ! atomic coordinates

! Local variables

   integer :: i 

   character(len=BUFLEN) :: infix   ! input operand tokens
   character(len=BUFLEN) :: postfix ! operand tokens in RPN form

   logical :: error

   ! check for null input:
   if( len_trim(maskstr) == 0 ) then
      mask(1:natom) = 0
      return
   end if

   ! put termination symbol ';' at the end of the string
   i = len_trim(maskstr) + 1
   maskstr(i:i) = ';'

   if (prnlev >= 1) then
      write(*,'("original : ==",A,"==")') maskstr(1:index(maskstr,';'))
   end if

   ! 1) preprocess input string expression
   call tokenize(maskstr, infix)
   if (prnlev >= 1) then
      write(*,'("tokenized: ==",A,"==")') infix(1:index(infix,';'))
   end if
   
   ! 2) construct postfix (RPN) notation
   call torpn(infix, postfix, error)
   if (error) call error1("atommask","unbalanced parantheses")
   if (prnlev >= 1) then
      write(*,'("postfix  : ==",A,"==")') postfix(1:index(postfix,';'))
   end if

   ! 3) evaluate postfix notation and return the result in mask array
   call eval(natom,nres,prnlev,igraph,isymbl,ipres,lbres,crd,postfix,mask)

end subroutine atommask

!*******************************************************************************
!
! Subroutine:  tokenize
!
! Description: 
!     
! Mark operand (anything enclosed between operators) tokens by enclosing them in
! [..]. This makes processing 'maskstr' to RPN easier. Operands are represented
! by atomic expressions such as [@CA], [:LYS,ARG], [:1-10@CB]
!              
!*******************************************************************************

subroutine tokenize(input, infix)
   implicit none
   
! Passed variables

   character(*)              :: input ! original mask string (last char is ;)
   character(*), intent(out) :: infix ! operand-tokenized string

! Local variables

   character(BUFLEN)  :: buffer   ! storage buffer
   character(1)       :: symbol   ! character buffer

   integer            :: i        ! position in 'infix'
   integer            :: p        ! position in 'input'
   integer            :: inplen   ! length of input string
   integer            :: j    
   integer            :: n

   logical            :: flag     ! general flag for logical control
   
   flag = .false.
   p = 1
   inplen = index(input,';')
   
   infix = ' '
   i = 1

   do p = 1, inplen

      ! skip blanks
      if (input(p:p) == ' ') cycle

      symbol = input(p:p)
      if (isOperand(symbol)) then
         if (.not.flag) then
            n = 1
            buffer(n:n) = '['
            n = n + 1
            flag = .true.
         end if
         buffer(n:n) = symbol
         n = n + 1
      else if (isOperator(symbol).or.(symbol=='(').or.  &
              (symbol==')').or.(symbol==';')) then
         if (flag) then
            buffer(n:n) = ']'
            n = n + 1
            if (buffer(2:2)==':'.and.index(buffer(1:n-1),'@')>0) then
               ! this expression has [:..@..] form and needs splitting
               infix(i:i) = '('
               i = i + 1
               do j=1,n-1
                  if (buffer(j:j) == '@') then
                     infix(i:i) = ']'
                     infix(i+1:i+1) = '&'
                     infix(i+2:i+2) = '['
                     i = i + 3
                  end if
                  infix(i:i) = buffer(j:j)
                  i = i + 1
               end do
               infix(i:i) = ')'
               i = i + 1
            else ! expression doesn't need spliting, just copy it out
               do j=1,n-1
                  infix(i:i) = buffer(j:j)
                  i = i + 1
               end do
            end if
            flag = .false.
         end if
         infix(i:i) = symbol
         i = i + 1
      else
         call error2("tokenize","illegal symbol in maskstr: ",symbol)
      end if
   end do
   
end subroutine tokenize

!*******************************************************************************
!
! Subroutine:  torpn
!
! Description: 
!
! Converts infix expression to RPN. ';' marks the end of infix. The basic stack
! operations are implemented by subprograms init_stack, empty_stack, push_stack,
! and pop_stack
!              
!*******************************************************************************

subroutine torpn(infix, postfix, error)
   implicit none

! Passed variables

   character(*) infix   ! input tokenized string
   character(*) postfix ! RPN-output tokenized string

   logical error        ! if error was found in infix

! Local variables

   character(1) symbol  ! character buffer
   character(1) topop   ! character buffer popped from stack

   integer, parameter :: maxstack=BUFLEN ! max stack size

   character(1) stack(maxstack)       ! stack character array

   integer :: top         ! top of stack location
   integer :: i           ! position in postfix
   integer :: p           ! position in infix
   integer :: inplen      ! length of infix string

   logical done

   ! Initialize an empty stack and blank our RPN
   call init_stack(top)
   postfix = ' '

   error = .false.
   i = 1
   p = 1
   inplen = index(infix,';') - 1
   
   ! while not end of 'infix' do:
   do p=1,inplen
      symbol = infix(p:p)
      if (isOperand(symbol).or.symbol=='['.or.symbol==']') then
         ! copy everything which is part of the operand 
         ! (icluding '[', ']') to 'postfix'
         postfix(i:i) = symbol
         i = i + 1
      else if (symbol.eq.'(') then
         ! left parentheses, push it onto stack
         call push_stack(stack, top, maxstack, symbol)
      else if (symbol.eq.')') then
         ! if symbol is right parentheses: pop symbols from stack until
         ! the corresponding left parentheses appears (adding these 
         ! symbols to the 'postfix' expression at possition 'i') or until
         ! stack becomes empty. 'error' is returned as true if no left 
         ! parentheses is found.
10       continue
            if (empty_stack(top)) then
               error = .true.
            else
               call pop_stack(stack, top, symbol)
               if (symbol.ne.'(') then
                  postfix(i:i) = symbol
                  i = i + 1
               end if
            end if
         if (.not.error.and.symbol.ne.'(') goto 10

      else if (isOperator(symbol)) then
         ! if symbol is an operator: pop operators from stack (max sized 
         ! maxstack, top is top of the stack) until the stack becomes empty
         ! or an operator appears on the top of the stack whose priority
         ! is less than or equal to that of operator symbol. symbol is then
         ! pushed onto the stack.
         !    topop : operator on top of stack
         !    done  : signals completion of stack popping

         done = .false.
         
         ! repeat the following until done popping
20       continue
            if (empty_stack(top)) then
               done = .true.
            else
               call pop_stack(stack,top,topop)
               if (priority(symbol).le.priority(topop)) then
                  postfix(i:i) = topop
                  i = i + 1
               else
                  call push_stack(stack, top, maxstack, topop)
                  done= .true.
               end if
            end if
         if (.not. done) goto 20
         
         call push_stack(stack,top,maxstack,symbol)

      else
         ! illegal character -- skip it
         write(*,'(A1,"is illegal character -- ignored")') symbol
      end if
   end do 
   ! end of 'infix' reached. If no error detected, pop any
   ! operands on the stack and put them in 'postfix'
   
   if (.not.error) then
30    if (.not. empty_stack(top)) then
         call pop_stack(stack, top, symbol)
         postfix(i:i) = symbol
         i = i + 1
         goto 30
      end if
   end if
   postfix(i:i) = ';'   ! terminate the 'postfix' string
   
contains

!*******************************************************************************
!
! Function: priority 
!
! Description: find priority of a given operator
!              
!*******************************************************************************

integer function priority(op)
   implicit none

   character(1) op  ! operator
   
   if (op.eq.'<'.or.op.eq.'>') then
      priority = 6
   else if (op.eq.'!') then
      priority = 5
   else if (op.eq.'&') then
      priority = 4
   else if (op.eq.'|') then
      priority = 3
   else if (op.eq.'(') then
      priority = 2
   else if (op.eq.';') then
      priority = 1
   else
      call error2("priority","unknown operator on stack: ",op)
   end if

end function priority

end subroutine torpn

!*******************************************************************************
!
! Subroutine:  eval
!
! Description: 
!
! Take RPN mask string, split into elementary selections, and translate into
! a mask integer array
!              
!*******************************************************************************

subroutine eval(natom, nres, prnlev, igraph, isymbl, ipres, lbres,  &
                crd, postfix, mask)
   implicit none
   
! Passed variables

   integer           :: natom      ! number of atoms
   integer           :: nres       ! number of residues
   integer           :: ipres(*)   ! residue pointer array
   integer           :: mask(*)    ! mask selection array
   integer           :: prnlev     ! print level

   character(len=4)  :: igraph(*)  ! atom names
   character(len=4)  :: isymbl(*)  ! atom types
   character(len=4)  :: lbres(*)   ! residue labels

   double precision  :: crd(*)       ! coordinate array

   character(*) postfix              ! RPN, tokenized mask operands

! Local variables

   integer, parameter :: maxstack=BUFLEN
   
   integer  :: i, j, p   ! string indexes
   integer  :: inplen    ! length of postfix
   integer  :: top       ! location of top of stack
   integer  :: astat     ! allocate status
   integer  :: nselatom  ! number of selected atoms

   character(len=BUFLEN) buffer, diststr
   character(1) symbol
   integer, dimension(:,:), allocatable :: stack
   integer, dimension(:), allocatable :: mask1, mask2
   
   ! allocate stack for (at most) 'maxstack' atom mask() arrays
   allocate(stack(maxstack,natom), stat=astat)
   if (astat/=0) call error1("eval","stack allocation error")
   
   call init_evalstack(top)  ! initialize Stack
   nselatom = 0
   i = 1
   p = 1
   inplen = index(postfix,';') - 1   ! up to (but not including) ';'
   
   do p=1,inplen
      symbol = postfix(p:p)
      if (symbol.eq.'[') then        ! 'operand' begins here
         i = 1
      else if (symbol.eq.']') then   ! 'operand' is completed
         buffer(i:i) = ';'           ! terminate elementary mask expression
         ! if operand is immediately followed by '<' or '>', it must be a 
         ! distance operand (: or @ followd by a number), and we should not 
         ! call selectElemMask() which is for mask operands, but store it 
         ! as it is and deal with it later in noneq()
         if (postfix(p+1:p+1).eq.'<'.or.postfix(p+1:p+1).eq.'>') then
            diststr = ''; diststr = buffer(1:i)  ! e.g. :3.1;  or @4.3;
            if (prnlev >= 3) then
               write(*,'("buffer: ",A)') trim(buffer)
               write(*,'("distance saved: ",A)') trim(diststr)
            end if
         else
            ! mask was allocated in a calling routine already (e.g. ambmask.f)
            ! print *, "buffer=",buffer(1:index(buffer,';'))
            call selectElemMask(natom,nres,igraph,isymbl,ipres,lbres,  &
                                buffer,mask)
            ! push the mask(1:natom) array to the stack
            call push_evalstack(stack,top,maxstack,natom,mask)
            if (prnlev >= 3) then
               write(*,'("pushed (top=",I2,")")') top
               write(*,'(50I1)') (mask(j),j=1,natom)
            end if
         end if
      else if (isOperand(symbol)) then   ! operand is a part inside [...]
         buffer(i:i) = symbol
         i = i + 1
      else if (symbol.eq.'&'.or.symbol.eq.'|') then
         allocate(mask1(natom), stat=astat)
         if (astat/=0) call error1("eval","mask1 (&,|) allocation error")
         allocate(mask2(natom), stat=astat)
         if (astat/=0) call error1("eval","mask2 (&,|) allocation error")
         call pop_evalstack(stack,top,maxstack,natom,mask1)
         call pop_evalstack(stack,top,maxstack,natom,mask2)
         if (prnlev >= 3) then
            write(*,'("popped (top=",2I2,")")') top+1, top
            write(*,'(50I1)') (mask2(j),j=1,natom)
            write(*,'(50I1)') (mask1(j),j=1,natom)
         end if
         ! binop performs the operation and returns the result in mask
         ! you can release both mask1 and mask2 upon return from binop()
         if (allocated(mask2).and.allocated(mask1)) then
            call binop(symbol, mask2, mask1, natom, mask)
            if (prnlev >= 3) then
               write(*,'("binary operation =",A1,"= (2 pops above)")') symbol
               write(*,'(50I1)') (mask(j),j=1,natom)
            end if
         else
            call error1("eval","illegal binary (&,|) operation")
         end if
         call push_evalstack(stack,top,maxstack,natom,mask)
         if (prnlev >= 3) then
            write(*,'("pushed (top=",I2,")")') top
            write(*,'(50I1)') (mask(j),j=1,natom)
         end if
         deallocate(mask1, stat=astat)
         if (astat/=0) call error1("eval","mask1 (&,|) deallocation error")
         deallocate(mask2, stat=astat)
         if (astat/=0) call error1("eval","mask2 (&,|) deallocation error")
      else if (symbol.eq.'<'.or.symbol.eq.'>') then
         allocate(mask2(natom), stat=astat)
         if (astat/=0) call error1("eval","mask2 (<,>) allocation error")
         ! the distance and if 'by (:)' or 'by (@)' is saved in 'diststr'
         ! get the mask from which we search 'distance' away
         call pop_evalstack(stack,top,maxstack,natom,mask2)
         if (prnlev >= 3) then
            write(*,'("popped (top=",I2,")")') top
            write(*,'(50I1)') (mask2(j),j=1,natom)
            write(*,'("distance restored ",A)') trim(diststr)
         end if
         ! noneq performs the operation and returns the result in mask
         ! you can release mask2 upon return from noneq()
         if (allocated(mask2)) then
            call noneq(symbol, mask2, diststr, natom, nres, ipres, crd, mask)
            if (prnlev >= 3) then
               write(*,'("binary operation =",A1,"= (2 pops above)")') symbol
               write(*,'(50I1)') (mask(j),j=1,natom)
            end if
         else
            call error1("eval","illegal binary (<,>) operation")
         end if
         call push_evalstack(stack,top,maxstack,natom,mask)
         if (prnlev >= 3) then
            write(*,'("pushed (top=",I2,")")') top
            write(*,'(50I1)') (mask(j),j=1,natom)
         end if
         deallocate(mask2, stat=astat)
         if (astat/=0) call error1("eval","mask2 (<,>) deallocation error")
      else if (symbol.eq.'!') then
         allocate(mask1(natom), stat=astat)
         if (astat/=0) call error1("eval","mask1 allocation error")
         call pop_evalstack(stack,top,maxstack,natom,mask1)
         if (prnlev >= 3) then
            write(*,'("popped (top=",I2,")")') top
            write(*,'(50I1)') (mask1(j),j=1,natom)
         end if
         if (allocated(mask1)) then
            call neg(mask1,natom,mask)
            if (prnlev >= 3) then
               write(*,'("unary operation =",A1,"= (1 pop above)")') symbol
               write(*,'(50I1)') (mask(j),j=1,natom)
            end if
         else
            call error1("eval","illegal unary neg operation")
         end if
         call push_evalstack(stack,top,maxstack,natom,mask)
         if (prnlev >= 3) then
            write(*,'("pushed (top=",I2,")")') top
            write(*,'(50I1)') (mask(j),j=1,natom)
         end if
         deallocate(mask1, stat=astat)
         if (astat/=0) call error1("eval","mask1 deallocation error")
      else
         call error1("eval","unknown symbol while evaluating RPN")
      end if
   end do

   ! Mask should be the resulting mask, and this should also free up 
   ! the last item on stack
   call pop_evalstack(stack,top,maxstack,natom,mask) 
   if (prnlev >= 2) then
      write(*,'("popped (top=",I2,")")') top
      write(*,'("this should be the resulting mask")')
      write(*,'(50I1)') (mask(j),j=1,natom)
   end if

   ! deallocate stack, we do not need it anymore
   deallocate(stack, stat=astat)
   if (astat/=0) call error1("eval","stack deallocation error")
   
   if (prnlev >= 1) then
      do j=1,natom
         if (mask(j).eq.1) nselatom = nselatom + 1
      end do
      if (nselatom == 0) then
         write(6,'("Warning: no atoms selected")')
      else 
         if (natom < 1000000) then
           write(6,'(I6," atoms selected")') nselatom
         else
           write(6,'(I9," atoms selected")') nselatom
         end if
      end if
   end if

end subroutine eval


!*******************************************************************************
!
! Subroutine:  selectElemMask
!
! Description: 
!
! Fill 'mask' array corresponding to elementary atom expression. Tedious
! analysis of elementary atom expressions is done here. It needs to find out
! whether it is a list of res numbers, res names, and likewise for atoms.
! It then parses these lists and converts them to a mask array
!
!*******************************************************************************

subroutine selectElemMask(natom, nres, igraph, isymbl, ipres, lbres, buffer, mask)
   implicit none

! Internal parameters 

   integer, parameter :: ALL = 0, NUMLIST = 1, NAMELIST = 2, ATYPELIST = 3

! Passed variables

   integer            :: natom     ! number of atoms
   integer            :: nres      ! number of residues
   integer            :: ipres(*)  ! residue pointers
   integer            :: mask(*)   ! mask array

   character (len=4)  :: igraph(*) ! atom names
   character (len=4)  :: isymbl(*) ! atom types
   character (len=4)  :: lbres(*)  ! res labels
   character(*)       :: buffer    ! buffer with elementary selections

! Local variables
   
   integer  :: i          
   integer  :: reslist
   integer  :: atomlist

   character(1) :: symbol ! character holder

   ! reset the mask for elementary selection (to nothing selected)
   do i=1,natom
      mask(i) = 0
   end do
   
   if (buffer(1:1).eq.':') then        ! residue mask expression
      if (buffer(2:2).eq.'*' ) then
         reslist = ALL
      else
         reslist = NUMLIST
         ! if there is a letter in it, it is NAMELIST
         do i=2,index(buffer,';')-1
            symbol = buffer(i:i)
            if (symbol.ge.'a'.and.symbol.le.'z'.or.  &
                symbol.ge.'A'.and.symbol.le.'Z') then
               reslist = NAMELIST
               exit
            end if
         end do
      end if
      select case (reslist)
         case (ALL)
            call all_select(mask, natom)
         case (NUMLIST)
            call residue_numlist(buffer(2:), mask, nres, ipres)
         case (NAMELIST) 
            call residue_namelist(buffer(2:), mask, nres, lbres, ipres)
      end select
   else if (buffer(1:1).eq.'@') then   ! atom mask expression
      ! because atom names can have digits, and even can start with a digit,
      ! we need to search the whole expression to decide whether it is
      ! an atom numlist or namelist and it is still ambiguous
      if (buffer(2:2).eq.'*') then
         atomlist = ALL
      else if (buffer(2:2).eq.'%') then
         atomlist = ATYPELIST      ! amber type list
      else
         atomlist = NUMLIST
         do i=2,index(buffer,';')-1
            symbol = buffer(i:i)
            if (symbol.ge.'a'.and.symbol.le.'z'.or.  &
                symbol.ge.'A'.and.symbol.le.'Z') then
               atomlist = NAMELIST
               exit
            end if
         end do
      end if
      select case (atomlist)
         case (ALL)
            call all_select(mask, natom)
         case (NUMLIST)
            call atom_numlist(buffer(2:), mask, natom)
         case (NAMELIST)
            call atom_namelist(buffer(2:), mask, natom, igraph)
         case (ATYPELIST)  ! there is '%' after '@'
            call atom_atypelist(buffer(3:), mask, natom, isymbl)
      end select
   else if (buffer(1:1).eq.'*') then   ! select all
      ! this is here just for compatibility with ptraj's notion of
      ! selecting all residues by '*' as opposed to ":*"
      call all_select(mask, natom)
   else
      call error1("selectElemMask",  &
           "elementary mask does not start with : or @")
   end if

end subroutine selectElemMask

!*******************************************************************************
!
! Subroutine:  binop
!
! Description: perform binary operation (& is and, | is or) on two mask arrays
!              
!*******************************************************************************

subroutine binop(op, mask2, mask1, natom, mask)
   implicit none
   
! Passed variables

   integer, intent(out) :: mask(*)  ! resulting mask
   integer, intent(in)  :: mask1(*) ! input masks to compare
   integer, intent(in)  :: mask2(*)
   integer, intent(in)  :: natom    ! number of atoms

   character(*)         :: op       ! operator

! Local variables

   integer  :: i
   
   ! reset resulting mask first:
   do i=1,natom
      mask(i) = 0
   end do
   
   select case (op)
      case ('&')
         do i=1,natom
            if (mask2(i).eq.1.and.mask1(i).eq.1) mask(i) = 1
         end do
      case ('|')
         do i=1,natom
            if (mask2(i).eq.1.or.mask1(i).eq.1) mask(i) = 1
         end do
      case default
         call error2("binop","unknown operator: ",op)
   end select
end subroutine binop

!*******************************************************************************
!
! Subroutine:  noneq
!
! Description: perform binary operation (<, >) on a mask array and distance
!              
!*******************************************************************************

subroutine noneq(op, mask2, diststr, natom, nres, ipres, crd, mask)
   implicit none

! Passed variables

   integer, intent(in)  :: mask2(*)       ! input mask
   integer, intent(in)  :: ipres(*)       ! residue pointers
   integer, intent(in)  :: natom          ! number of atoms
   integer, intent(in)  :: nres           ! number of residues

   integer, intent(out) :: mask(*)        ! resulting mask
   
   character(*) :: op                     ! operator
   character(*) :: diststr                ! distance string

   double precision, intent(in) :: crd(*) ! coordinate array

! Local variables

   integer :: i, j, k, ii, jj, ios

   double precision :: dist, d2, r2

   ! reset resulting mask first:
   do i=1,natom
      mask(i) = 0
   end do

   i = index(diststr,';') - 1
   if (i > 0) then
      read(diststr(2:i),*,iostat=ios) dist
      ! debug: write(*,'("dist=",F8.2)') dist
      if (ios > 0) call error1("noneq", "error parsing distance cutoff value")
      d2 = dist*dist
   else
      call error1("noneq", "error parsing distance cutoff value")
   end if

   ! First, find atoms closer than 'dist' (regardless of whether the
   ! selection is atom or residue based and regardless of the distance
   ! operator ('<' or '>') used).
   if (diststr(1:1).eq.'@'.or.diststr(1:1).eq.':') then
      do i=1,natom
         if (mask2(i) == 1) then
            ii = 3*i - 3
            do j=1,natom
               jj = 3*j - 3
               r2 = (crd(ii+1)-crd(jj+1))**2 +  &
                    (crd(ii+2)-crd(jj+2))**2 +  &
                    (crd(ii+3)-crd(jj+3))**2
               if (r2 < d2) mask(j) = 1
            end do ! j
         end if
      end do ! i
   end if

   if (diststr(1:1).eq.':') then    ! residue based cutoff
      ! if at least one atom in a residue is selected (from
      ! the previous step), select the whole residue
      do j = 1, nres   ! residue index
         do i = ipres(j), ipres(j+1)-1   ! atom index
            if (mask(i) == 1) then
               do k = ipres(j), ipres(j+1)-1
                  mask(k) = 1
               end do
               exit   ! exit i-th residue loop
            end if
         end do
      end do
   end if
   
   if (op.eq.'<') then
      continue  ! you're done
   else if (op.eq.'>') then
      ! Invert the selection: the best way of thinking about it is that
      ! if you want atoms/residues further away than 'dist', select
      ! atoms/residues that are closer than 'dist' and neg this selection.
      do i=1,natom
         if (mask(i) == 1) then
            mask(i) = 0
         else
            mask(i) = 1
         end if
      end do ! i
   else
      call error1("noneq","operator is not < or >")
   end if
   
   if (diststr(1:1).ne.'@'.and.diststr(1:1).ne.':') then
      call error1("noneq", "residue (<:,>:) or atom (<@,>@) " //  &
                  "based distance cutoff must be specified")
   end if

end subroutine noneq
   
!*******************************************************************************
!
! Subroutine:  neg
!
! Description: perform unary operation ! (negate) on a mask array
!              
!*******************************************************************************

subroutine neg(mask1, natom, mask)
   implicit none

! Passed variables

   integer, intent(in)  :: mask1(*)  ! input mask
   integer, intent(in)  :: natom     ! number of atoms
   integer, intent(out) :: mask(*)   ! output mask

! Local variables

   integer :: i

   do i=1,natom
      if (mask1(i).eq.1) then
         mask(i) = 0
      else
         mask(i) = 1
      end if
   end do
end subroutine neg


!********************************************************************************
!
! The routines below deal with parsing elementary expressions,
! which were obtained by the routine eval(..,postfix,mask), etc.
!
!********************************************************************************


!*******************************************************************************
!
! Subroutine:  residue_numlist
!
! Description: split residue number list to individual residue numbers
!              
!*******************************************************************************

subroutine residue_numlist(numlist, mask, nres, ipres)
   implicit none
   
! Passed variables

   character(*)  :: numlist   ! number list [:1,5,7], [:1-4], etc.

   integer       :: nres      ! number of residues
   integer       :: ipres(*)  ! residue pointers
   integer       :: mask(*)   ! output mask

! Local variables

   character(BUFLEN) :: buffer   ! buffer
   character(1)      :: symbol   ! character holder
   integer p, i, inplen, res1, res2, ios
   logical dash
   
   i = 1
   res1 = 1; res2 = 1
   dash = .false.
   inplen = index(numlist,';')-1
   do p=1,inplen
      symbol = numlist(p:p)
      if (symbol.ge.'0'.and.symbol.le.'9') then
         buffer(i:i) = symbol
         i = i + 1
      end if
      if (symbol.eq.','.or. p == inplen) then
         if (.not. dash) then
            !print *, "res1: >>",buffer(1:i-1),"<<"
            read(buffer(1:i-1),*,iostat=ios) res1
            ! check for read errors:
            if (ios > 0) call error1("residue_numlist",  &
                              "error reading residue number")
            ! check for res1 value not necessary: resnum_select()
            ! will not cause "out of array boundary"
            ! if (res1 < 1) res1 = 1
            ! if (res1 > nres) res1 = nres
            call resnum_select(res1, res1, mask, nres, ipres)
         else
            !print *, "res2: >>",buffer(1:i-1),"<<"
            read(buffer(1:i-1),*,iostat=ios) res2
            ! check for read errors:
            if (ios > 0) call error1("residue_numlist",  &
                              "error reading residue number")
            call resnum_select(res1, res2, mask, nres, ipres)
            dash = .false.
         end if
         i = 1
      else if (symbol.eq.'-') then
         !print *, "res1: >>",buffer(1:i-1),"<<"
         read(buffer(1:i-1),*,iostat=ios) res1
         ! check for read errors:
         if (ios > 0) call error1("residue_numlist",  &
                           "error reading residue number")
         dash = .true.
         i = 1
      end if
      if (.not. (symbol.ge.'0'.and.symbol.le.'9'.or.  &
                 symbol.eq.','.or.symbol.eq.'-') ) then
         call error2("residue_numlist","unknown symbol: ",symbol)
      end if
   end do

end subroutine residue_numlist

!********************************************************************************
!
! Subroutine:  residue_namelist
!
! Description: Split residue name list into individual residue names
!
!********************************************************************************

subroutine residue_namelist(namelist, mask, nres, lbres, ipres) 
   implicit none

! Passed variables

   character(*)  :: namelist    ! residue name list, i.e. [:ASP,ALA]
   character(4)  :: lbres(*)    ! residue labels

   integer       :: nres        ! number of residues
   integer       :: ipres(*)    ! residue pointers
   integer       :: mask(*)     ! output mask

! Local variables

   character(4) resname
   character(1) symbol
   integer p, i, inplen

   i = 1
   resname = "    "
   inplen = index(namelist,';')-1

   do p=1,inplen
      symbol = namelist(p:p)
      ! residue names may have: letters,digits, + and - (for ions),
      ! and '=' wildcard
      if (symbol.ge.'a'.and.symbol.le.'z'.or.  &
          symbol.ge.'A'.and.symbol.le.'Z'.or.  &
          symbol.ge.'0'.and.symbol.le.'9'.or.  &
          symbol.eq.'+' .or.symbol.eq.'-'.or.symbol.eq.'=') then
         if (i > 4) call error1("residue_namelist",  &
                    "residue name should not have more than 4 chars")
         resname(i:i) = symbol
         i = i + 1
      end if
      if (symbol.eq.','.or. p == inplen) then
         call resname_select(resname,mask,nres,lbres,ipres)
         i = 1
         resname = "    "  ! has to be cleared for next use
      end if
      if (.not.(symbol.ge.'a'.and.symbol.le.'z'.or.  &
                symbol.ge.'A'.and.symbol.le.'Z'.or.  &
                symbol.ge.'0'.and.symbol.le.'9'.or.  &
                symbol.eq.'+'.or.symbol.eq.'-'.or.   &
                symbol.eq.','.or.symbol.eq.'=') ) then
         call error2("residue_namelist","unknown symbol: ",symbol)
      end if
   end do

end subroutine residue_namelist

!********************************************************************************
!
! Subroutine:  atom_numlist
!
! Description: split atom numbers list into individual atom numbers
!
!********************************************************************************

subroutine atom_numlist(numlist, mask, natom)
   implicit none
   
! Passed variables

   character(*)  :: numlist ! number list i.e. [@1,3-6]
   integer       :: natom   ! number of atoms
   integer       :: mask(*) ! output mask

! Local variables

   character(len=BUFLEN) buffer
   character(1)  symbol
   integer p, i, inplen, at1, at2, ios
   logical dash

   i = 1
   at1 = 1
   at2 = 1
   dash = .false.
   inplen = index(numlist,';')-1
   
   do p=1,inplen
      symbol = numlist(p:p)
      if (symbol.ge.'0'.and.symbol.le.'9') then
         buffer(i:i) = symbol
         i = i + 1
      end if
      if (symbol.eq.','.or. p == inplen) then
         if (.not. dash) then
            read(buffer(1:i-1),*,iostat=ios) at1
            ! check for read errors:
            if (ios > 0) call error1("atom_numlist",  &
                              "error reading atom number")
            ! check at1 value not necessary, atnum_select()
            ! will not cause "out of array boundary"
            ! if (at1 < 1) at1 = 1
            ! if (at1 > natom) at1 = natom
            call atnum_select(at1, at1, mask, natom)
         else
            read(buffer(1:i-1),*,iostat=ios) at2
            ! check for read errors:
            if (ios > 0) call error1("atom_numlist",  &
                              "error reading atom number")
            call atnum_select(at1, at2, mask, natom)
            dash = .false.
         end if
         i = 1
      else if (symbol.eq.'-') then
         read(buffer(1:i-1),*,iostat=ios) at1
         ! check for read errors:
         if (ios > 0) call error1("atom_numlist",  &
                           "error reading atom number")
         dash = .true.
         i = 1
      end if
      if (.not. (symbol.ge.'0'.and.symbol.le.'9'.or.  &
                 symbol.eq.','.or.symbol.eq.'-') ) then
         call error2("atom_numlist","unknown symbol: ",symbol)
      end if
   end do

end subroutine atom_numlist

!********************************************************************************
!
! Subroutine:  atom_namelist
!
! Description: split atom name list into individual atom names
!
!********************************************************************************

subroutine atom_namelist(namelist, mask, natom, igraph)
! examples of atom namelists: [@CA,C,N] [@C5',C3',Na+,Cl-]
   implicit none

! Passed variables

   character(*)     :: namelist
   character(len=4) :: igraph(*)

   integer          :: natom
   integer          :: mask(*)

! Local variables

   character(4) atname
   character(1) symbol
   integer p, i, inplen
   
   i = 1
   atname = "    "
   inplen = index(namelist,';')-1

   do p=1,inplen
      symbol = namelist(p:p)
      ! atom names may have: letters, digits, primes (C5'),
      ! + and - (for ions Na+, Cl-), and '=' wildcard
      if (symbol.ge.'a'.and.symbol.le.'z'.or.  &
          symbol.ge.'A'.and.symbol.le.'Z'.or.  &
          symbol.ge.'0'.and.symbol.le.'9'.or.  &
          symbol.eq.'+' .or.symbol.eq.'-'.or.  &
          symbol.eq.''''.or.symbol.eq.'=') then
         if (i > 4) call error1("atom_namelist",  &
                    "atom name should not have more than 4 chars")
         atname(i:i) = symbol
         i = i + 1
      end if
      if (symbol.eq.','.or. p == inplen) then
         call atname_select(atname,mask,natom,igraph)
         i = 1
         atname = "    "  ! has to be cleared for next use
      end if
      if (.not.(symbol.ge.'a'.and.symbol.le.'z'.or.  &
                symbol.ge.'A'.and.symbol.le.'Z'.or.  &
                symbol.ge.'0'.and.symbol.le.'9'.or.  &
                symbol.eq.'+' .or.symbol.eq.'-'.or.  &
                symbol.eq.''''.or.symbol.eq.','.or.symbol.eq.'=') ) then
         call error2("atom_namelist","unknown symbol: ",symbol)
      end if
   end do

end subroutine atom_namelist

!********************************************************************************
!
! Subroutine:  atom_atypelist
!
! Description: Split amber atom type list into individual atom types
!
!********************************************************************************

subroutine atom_atypelist(atypelist, mask, natom, isymbl)
   implicit none

! Passed variables

   character(*) :: atypelist  ! atom type list. i.e. [@%C,CT]
   character(4) :: isymbl(*)

   integer      :: natom
   integer      :: mask(*)

! Local variables

   character(4) atatype
   character(1) symbol

   integer p, i, inplen
   
   i = 1
   atatype = "    "   ! blank ATtom Amber TYPE 'atatype'
   inplen = index(atypelist,';')-1

   do p=1,inplen
      symbol = atypelist(p:p)
      ! atomic element types can only be letters, digits and asterisk [*]
      if (symbol.ge.'a'.and.symbol.le.'z'.or.  &
          symbol.ge.'A'.and.symbol.le.'Z'.or.  &
          symbol.ge.'0'.and.symbol.le.'9'.or.  &
          symbol.eq.'*') then
         if (i > 4) call error1("atom_atypelist",  &
                    "amber type should not have more than 4 chars")
         atatype(i:i) = symbol
         i = i + 1
      end if
      if (symbol.eq.','.or. p == inplen) then
         call atatype_select(atatype,mask,natom,isymbl)
         i = 1
         atatype = "    "  ! has to be cleared for next use
      end if
      if (.not.(symbol.ge.'a'.and.symbol.le.'z'.or.  &
                symbol.ge.'A'.and.symbol.le.'Z'.or.  &
                symbol.ge.'0'.and.symbol.le.'9'.or.  &
                symbol.eq.'*'.or. symbol.eq.',') ) then
         call error2("atom_atypelist","unknown symbol: ",symbol)
      end if
   end do

end subroutine atom_atypelist

!********************************************************************************
!
! Subroutine:  resnum_select
!
! Description: set mask array to 1 for the number range (res1-res2) of residues
!
!********************************************************************************

subroutine resnum_select(res1, res2, mask, nres, ipres) 
   implicit none

! Passed variables
   integer  :: res1
   integer  :: res2
   integer  :: nres
   integer  :: mask(*)
   integer  :: ipres(*)

! Local variables

   integer i, j

   do i=1,nres
      if (i >= res1 .and. i <= res2) then
         do j=ipres(i),ipres(i+1)-1
            mask(j) = 1
         end do
      end if
   end do
end subroutine resnum_select

!********************************************************************************
!
! Subroutine:  resname_select
!
! Description: set mask array to 1 for a residue with name 'resname'
!
!********************************************************************************

subroutine resname_select(resname, mask, nres, lbres, ipres)
   implicit none

! Passed variables

   character (len=4) :: resname
   character (len=4) :: lbres(*)

   integer :: mask(*)
   integer :: nres
   integer :: ipres(*)

! Local variables
   integer i, j

   do i=1,nres
      if (isNameMatch(lbres(i), resname)) then
         do j=ipres(i),ipres(i+1)-1
            mask(j) = 1
         end do
      end if
   end do
end subroutine resname_select

!********************************************************************************
!
! Subroutine:  all_select
!
! Description: Select all atoms
!
!********************************************************************************

subroutine all_select(mask, natom)
   implicit none
   
! Passed variables

   integer mask(*)
   integer natom

! Local variables

   integer j

   do j=1,natom
      mask(j) = 1
   end do
end subroutine all_select

!********************************************************************************
!
! Subroutine:  atnum_select
!
! Description: Select atoms in atom number range at1 - at2
!
!********************************************************************************

subroutine atnum_select(at1, at2, mask, natom)
   implicit none

! Passed variables

   integer at1
   integer at2
   integer natom
   integer mask(*)

! Local variables

   integer j
   
   do j=1,natom
      if (j >= at1 .and. j <= at2) mask(j) = 1
   end do
end subroutine atnum_select

!********************************************************************************
!
! Subroutine:  atname_select
!
! Description: Select atoms with atom name atname
!
!********************************************************************************

subroutine atname_select(atname, mask, natom, igraph) 
   implicit none

! Passed variables

   character (len=4) atname 
   character (len=4) igraph(*)

   integer mask(*)
   integer natom 

! Local variables

   integer j

   do j=1,natom
      if (isNameMatch(igraph(j), atname)) mask(j) = 1
   end do
end subroutine atname_select

!********************************************************************************
!
! Subroutine:  atatype_select
!
! Description: Select atoms with atom type 'atatype'
!
!********************************************************************************

subroutine atatype_select(atatype,mask,natom,isymbl)
   implicit none

! Passed variables

   character (len=4) atatype 
   character (len=4) isymbl(*)
   integer mask(*)
   integer natom
   
! Local variables

   integer j

   do j=1,natom
      if (isNameMatch(isymbl(j), atatype)) mask(j) = 1
   end do
end subroutine atatype_select

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Stack processing subprograms
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! There are two sets of stack routines:
! 1) First ones (init_stack, empty_stack, push_stack, pop_stack) are 
!    used by torpn() for conversion of original (infix) string into 
!    RPN (postfix). Stack is dimensioned statically (in routine torpn)
!    and can only hold 1 character in stack entry.
! 2) Second set (init_evalstack, push_evalstack, pop_evalstack)
!    is used in evaluating RPN (postfix) expression. Its size 
!    is also dimensioned statically (in routine eval) but each
!    stack entry holds 'mask' integer array, which is 'natom' ints
!    in size (that size is dimensioned dynamically in routine
!    eval).

!********************************************************************************
!
! Subroutine:  init_stack
!
! Description: initialize an empty stack of character strings
!
!********************************************************************************

subroutine init_stack(top)
   implicit none
   integer top
   
   top = 0

end subroutine init_stack

!********************************************************************************
!
! Function:  empty_stack
!
! Description: check if a set of character strings is empty
!
!********************************************************************************

function empty_stack(top)
   implicit none
   logical empty_stack
   integer top
   
   empty_stack = (top.eq.0)
   
end function empty_stack

!********************************************************************************
!
! Subroutine:  push_stack
!
! Description: push character symbol on the stack
!
!********************************************************************************

subroutine push_stack(stack, top, maxstack, symbol)
   implicit none
   integer top, maxstack
   character(1) stack(maxstack), symbol
   
   if (top.eq.maxstack) then
      call error1("push_stack","stack overflow:" //  &
                  "increase parameter maxstack in findmask.f::torpn()")
   else
      top = top + 1
      stack(top) = symbol
   end if
end subroutine push_stack

!********************************************************************************
!
! Subroutine:  pop_stack
!
! Description: Pops the top symbol from stack
!
!********************************************************************************

subroutine pop_stack(stack, top, symbol)
   implicit none
   integer top
   character(1) stack(*), symbol
   
   if (top.eq.0) then
      call error1("pop_stack","stack underflow")
   else 
      symbol = stack(top) 
      top = top - 1
   end if
end subroutine pop_stack


!********************************************************************************
!
! Subroutine:  init_evalstack
!
! Description: Initializes an empty stack (by setting top to zero)
!
!********************************************************************************

subroutine init_evalstack(top)
   implicit none
   integer top
   
   top = 0

end subroutine init_evalstack

!********************************************************************************
!
! Subroutine:  push_evalstack
!
! Description: pushes mask array on the stack
!
!********************************************************************************

subroutine push_evalstack(stack, top, maxstack, natom, mask)
   implicit none
   integer top, maxstack, natom, j
   integer, dimension(maxstack,natom) :: stack
   integer, dimension(natom) :: mask
   
   if (top.eq.maxstack) then
      call error1("push_evalstack","stack overflow:" //   &
                  "increase parameter maxstack in findmask.f::eval()")
   else
      top = top + 1
      do j=1,natom
         stack(top,j) = mask(j)
      end do
   end if
end subroutine push_evalstack

!********************************************************************************
!
! Subroutine:  pop_evalstack
!
! Description: pops top mask array from stack
!
!********************************************************************************

subroutine pop_evalstack(stack, top, maxstack, natom, mask)
   implicit none
   integer top, maxstack, natom, j
   integer, dimension(maxstack,natom) :: stack
   integer, dimension(natom) :: mask
   
   if (top.eq.0) then
      call error1("pop_evalstack","stack underflow")
   else
      do j=1,natom
         mask(j) = stack(top,j)
      end do
      top = top - 1
   end if
end subroutine pop_evalstack


!********************************************************************************
!
! Subroutine:  error1
!
! Description: Report error and quit program; takes 2 string arguments
!
!********************************************************************************

subroutine error1(subname, message)
   use pmemd_lib_mod, only : mexit

   character(*) subname
   character(*) message
   
   ! report which module (filename) and which routine error occurred in
   write(*,'("Error in group input::atommask.f::",A)') trim(subname)
   ! be more specific about the error
   write(*,'(A)') trim(message)
   ! quit the program gracefully
   call mexit(6,1)
end subroutine error1

!********************************************************************************
!
! Subroutine:  error2
!
! Description: report error and quit the program; takes 3 string arguments
!
!********************************************************************************

subroutine error2(subname, message1, message2)
   use pmemd_lib_mod, only : mexit

   character(*) subname
   character(*) message1
   character(*) message2
   
   ! report which module (filename) and which routine error occurred in
   write(*,'("Error in group input::atommask.f::",A)') trim(subname)
   ! be more specific about the error
   write(*,'(2A)') trim(message1), trim(message2)
   ! quit the program gracefully
   call mexit(6,1)
end subroutine error2

!********************************************************************************
!
! Function: isOperand 
!
! Description: 
!
! Logical-valued function that determines if symbol is an operand:
! operand can be: 
! [a-zA-Z0-9] is for names and numbers, ':' is residue, '@' is atom,
! '*' is everything or part of amber type, '=' is for a wildcard, 
! '-' for atom/residue range or part of ion name, ''' for some atom
! names (in nucleic acids), ',' for atom number enumeration, 
! '.' is for decimal point, '%' is for amber atom type, 
! '+' is for positive ions
!
!********************************************************************************

function isOperand(symbol)
   implicit none
   character(*) symbol
   logical isOperand
   
   isOperand = .false.
   if ((symbol.ge.'A').and.(symbol.le.'Z').or.  &
       (symbol.ge.'a').and.(symbol.le.'z').or.  &
       (symbol.ge.'0').and.(symbol.le.'9').or.  &
       (symbol.eq.':').or.(symbol.eq.'@').or.   &
       (symbol.eq.'*').or.(symbol.eq.'=').or.   &
       (symbol.eq.'-').or.(symbol.eq.',').or.   &
       (symbol.eq.'.').or.(symbol.eq.'%').or.   &
       (symbol.eq.'+').or.(symbol.eq.'''')) isOperand = .true.
       
end function isOperand

!*******************************************************************************
!
! Function:  isOperator
!
! Description: Logical-valued function that determines if symbol is an operator
!
!*******************************************************************************

function isOperator(symbol)
   implicit none
   character(*) symbol
   logical isOperator
   
   isOperator = .false.
   if (symbol.eq.'<') isOperator = .true.
   if (symbol.eq.'>') isOperator = .true.
   if (symbol.eq.'!') isOperator = .true.
   if (symbol.eq.'&') isOperator = .true.
   if (symbol.eq.'|') isOperator = .true.

end function isOperator

!*******************************************************************************
!
! Function:  isNameMatch
!
! Description: compare two strings representing atom or amber atom names
!
!*******************************************************************************

function isNameMatch(s1,s2)
   implicit none
   logical isNameMatch
   character(4) s1, s2    ! s1 is from prmtop, s2 is from maskstr
   character(4) s1u, s2u  ! left adjusted uppercase strings s1 and s2
   integer i, idx         ! indices used for '=' wildcard
   
   ! Both strings should be already lef justified but we
   ! do it anyway to be sure (the comparison depends on it).
   ! Also, we create new (left adjusted and uppercased) strings,
   ! such that original string arrays (e.g. igraph) do not get
   ! modified by the changes we make here
   s1u = " "; s2u = " "   ! clear strings
   call str2upper(adjustl(s1),s1u)
   call str2upper(adjustl(s2),s2u)
   
   ! deal with '=' wildcard
   idx = index(s2u,'=')
   if (idx > 0 .and. idx <5) then
      do i = idx, len(s2)
         s2u(i:i) = ' '   ! blank all symbols following '=' (including '=')
         s1u(i:i) = ' '
      end do
   end if
   
   if (s1u.eq.s2u) then
      isNameMatch = .true.
   else
      isNameMatch = .false.
   end if

end function isNameMatch

!********************************************************************************
!
! Subroutine:  str2upper
!
! Description: converts string to upper-case
!
!********************************************************************************

subroutine str2upper(str,upper)
   implicit none
   character(len=*), intent(in) :: str
   character(len=*), intent(out) :: upper
   integer :: i
   
   ! not sure if this will work on all kinds of machines
   ! a more general way to convert to uppercase may be needed
   integer, parameter :: offset = ichar('A') - ichar('a')
   do i=1,len_trim(str)
      if (lge(str(i:i),'a') .and. lle(str(i:i),'z')) then
         upper(i:i) =  char(ichar(str(i:i)) + offset)
      else
         upper(i:i) = str(i:i)
      end if
   end do
end subroutine str2upper

end module findmask_mod
