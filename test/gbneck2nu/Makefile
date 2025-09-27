test:
	-cd 1hji && ./Run.1hji_md
	-cd modified_nu/3p4a/ && ./Run.3p4a
	-export TESTsander='../pmemd'; cd 1hji && ./Run.1hji_md
	-export TESTsander='../../pmemd'; cd modified_nu/3p4a/ && ./Run.3p4a 

test.sander:
	-cd 1hji && ./Run.1hji_md
	-cd modified_nu/3p4a/ && ./Run.3p4a

test.pmemd:
	-export TESTsander='../pmemd'; cd 1hji && ./Run.1hji_md
	-export TESTsander='../../pmemd'; cd modified_nu/3p4a/ && ./Run.3p4a 

test.cuda:
	-export TESTsander='../pmemd.cuda'; cd 1hji && ./Run.1hji_md
	-export TESTsander='../../pmemd.cuda'; cd modified_nu/3p4a/ && ./Run.3p4a 
