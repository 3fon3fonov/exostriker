cd /home/beust/swiftori/hjs/tools
gen_ae_hjs <<!
1   !   AU & Years
plhjs_5.in
654876543
5 !  # of first set of test particles
0 -1 0 0 0 ! orbct : tp's orbit bodies 2
0. 0.1  ! Ecc's
3.  ! Max inc
1.5 3.5 ! Amin & Amax
0.     ! Power index  
7 !  # of second set of test particles
0 0 0 -1 0 ! orbct : tp's orbit bodies 4
0. 0.1  ! Ecc's
3.  ! Max inc
1.5 3.5 ! Amin & Amax
0.     ! Power index 
10 !  # of third set of test particles
-1 -1 0 0 0 ! orbct : tp's orbit bodies 1+2
0. 0.1  ! Ecc's
3.  ! Max inc
1.5 3.5 ! Amin & Amax
0.     ! Power index 
20 !  # of fourth set of test particles
-1 -1 -1 0 0 ! orbct : tp's orbit bodies 1+2+3
0. 0.1  ! Ecc's
3.  ! Max inc
1.5 3.5 ! Amin & Amax
0.     ! Power index 
10 !  # of fifth set of test particles
0 0 0 -1 -1 ! orbct : tp's orbit bodies 4+5
0. 0.1  ! Ecc's
3.  ! Max inc
1.5 3.5 ! Amin & Amax
0.     ! Power index 
7 !  # of sixth set of test particles
0 0 -1 0 0 ! orbct : tp's orbit bodies 3
0. 0.1  ! Ecc's
3.  ! Max inc
1.5 3.5 ! Amin & Amax
0.     ! Power index 
30 !  # of seventh set of test particles
-1 -1 -1 -1 -1 ! orbct : tp's orbit all bodies
0. 0.1  ! Ecc's
3.  ! Max inc
1.5 3.5 ! Amin & Amax
0.     ! Power index
0     ! 8th set of tp's => dummy
tphjs_5.in
!
exit
