$begin 'Profile'
	$begin 'ProfileGroup'
		MajorVer=2025
		MinorVer=2
		Name='Solution Process'
		$begin 'StartInfo'
			I(1, 'Start Time', '11/11/2026 21:50:39')
			I(1, 'Host', 'AEM-LAB1')
			I(1, 'Processor', '48')
			I(1, 'OS', 'NT 10.0')
			I(1, 'Product', 'HFSS Version 2025.2.0')
		$end 'StartInfo'
		$begin 'TotalInfo'
		$end 'TotalInfo'
		GroupOptions=8
		TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 1, \'Executing From\', \'C:\\\\Program Files\\\\ANSYS Inc\\\\v252\\\\AnsysEM\\\\HFSSCOMENGINE.exe\')', false, true)
		$begin 'ProfileGroup'
			MajorVer=2025
			MinorVer=2
			Name='HPC'
			$begin 'StartInfo'
				I(1, 'Type', 'Manual')
				I(1, 'Distribution Types', 'Variations, Frequencies, Domain Solver, Transient Excitations, Direct Solver, Iterative Solver, Solver, Mesh Assembly, Mesher')
				I(1, 'MPI Vendor', 'Intel')
				I(1, 'MPI Version', '2021')
			$end 'StartInfo'
			$begin 'TotalInfo'
				I(0, ' ')
			$end 'TotalInfo'
			GroupOptions=0
			TaskDataOptions(Memory=8)
			ProfileItem('Two Level', 0, 0, 0, 0, 0, 'I(1, 0, \'Disabled\')', false, true)
			ProfileItem('Machine', 0, 0, 0, 0, 0, 'I(6, 1, \'Name\', \'AEM-LAB1\', 1, \'Memory\', \'256 GB\', 3, \'RAM Limit\', 90, \'%f%%\', 2, \'Tasks\', 1, false, 2, \'Cores\', 5, false, 1, \'Free Disk Space\', \'310 GB\')', false, true)
		$end 'ProfileGroup'
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 1, \'Allow off core\', \'True\')', false, true)
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 1, \'Solution Basis Order\', \'1\')', false, true)
		ProfileItem('Design Validation', 0, 0, 0, 0, 0, 'I(1, 0, \'Elapsed time : 00:00:00 , HFSS ComEngine Memory : 231 M\')', false, true)
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'Perform full validations with standard port validations\')', false, true)
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
		$begin 'ProfileGroup'
			MajorVer=2025
			MinorVer=2
			Name='Initial Meshing'
			$begin 'StartInfo'
				I(1, 'Time', '11/11/2026 21:50:41')
			$end 'StartInfo'
			$begin 'TotalInfo'
				I(1, 'Elapsed Time', '00:18:15')
			$end 'TotalInfo'
			GroupOptions=4
			TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
			ProfileItem('Mesh', 608, 0, 781, 0, 5260000, 'I(3, 2, \'Tetrahedra\', 4183852, false, 1, \'Type\', \'TAU\', 2, \'Cores\', 5, false)', true, true)
			ProfileItem('Coarsen', 341, 0, 341, 0, 5260000, 'I(1, 2, \'Tetrahedra\', 1153984, false)', true, true)
			ProfileItem('Lambda Refine', 26, 0, 29, 0, 1128132, 'I(1, 2, \'Tetrahedra\', 1175399, false)', true, true)
			ProfileItem('Simulation Setup', 22, 0, 22, 0, 2765172, 'I(2, 2, \'Tetrahedra\', 774586, false, 1, \'Disk\', \'0 Bytes\')', true, true)
			ProfileItem('Port Adapt', 14, 0, 13, 0, 2779740, 'I(2, 2, \'Tetrahedra\', 774586, false, 1, \'Disk\', \'76.6 KB\')', true, true)
			ProfileItem('Port Refine', 71, 0, 77, 0, 1467036, 'I(1, 2, \'Tetrahedra\', 1170663, false)', true, true)
		$end 'ProfileGroup'
		$begin 'ProfileGroup'
			MajorVer=2025
			MinorVer=2
			Name='Adaptive Meshing'
			$begin 'StartInfo'
				I(1, 'Time', '11/11/2026 22:08:57')
			$end 'StartInfo'
			$begin 'TotalInfo'
			$end 'TotalInfo'
			GroupOptions=4
			TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
			$begin 'ProfileGroup'
				MajorVer=2025
				MinorVer=2
				Name='Adaptive Pass 1'
				$begin 'StartInfo'
					I(1, 'Frequency', '5.8GHz')
				$end 'StartInfo'
				$begin 'TotalInfo'
					I(0, ' ')
				$end 'TotalInfo'
				GroupOptions=0
				TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
				ProfileItem(' ', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
				ProfileItem('Simulation Setup ', 40, 0, 41, 0, 2756784, 'I(2, 2, \'Tetrahedra\', 774732, false, 1, \'Disk\', \'185 KB\')', true, true)
				ProfileItem('Matrix Assembly', 98, 0, 176, 0, 7584348, 'I(3, 2, \'Tetrahedra\', 774732, false, 2, \'Lumped ports\', 1, false, 1, \'Disk (output)\', \'301 MB\')', true, true)
				ProfileItem('Matrix Solve', 172, 0, 557, 0, 16464240, 'I(5, 1, \'Type\', \'MCS\', 2, \'Cores\', 5, false, 2, \'Matrix size\', 4044339, false, 3, \'Matrix bandwidth\', 17.4797, \'%5.1f\', 1, \'Disk\', \'15.4 MB\')', true, true)
				ProfileItem('Field Recovery', 39, 0, 66, 0, 16464240, 'I(2, 2, \'Excitations\', 1, false, 1, \'Disk\', \'97.6 MB\')', true, true)
				ProfileItem('Data Transfer', 24, 0, 6, 0, 250296, 'I(1, 0, \'Adaptive Pass 1\')', true, true)
				$begin 'ProfileGroup'
					MajorVer=2025
					MinorVer=2
					Name='APIPms'
					$begin 'StartInfo'
						I(1, 'Timesinceepock', '1794415452')
					$end 'StartInfo'
					$begin 'TotalInfo'
						I(0, ' ')
					$end 'TotalInfo'
					GroupOptions=16
					TaskDataOptions(Memory=8)
					ProfileItem('solverinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Solvertype\', \'shared_memory\', 1, \'Precision\', \'mixed\', 1, \'Solversymmetry\', \'complex_sym\', 1, \'Matrixdim\', \'4044339\', 1, \'Matrixbw\', \'17.481300\', 1, \'Matrixnnz\', \'70700202\', 1, \'Rootdim\', \'5512\', 1, \'Mathtype\', \'mkl\', 1, \'Mpitasks\', \'1\', 1, \'Threadspertasks\', \'0\')', false, true)
					ProfileItem('sysinfo', 0, 0, 0, 0, 0, 'I(12, 1, \'Os\', \'win\', 1, \'Cpuid\', \'Intel(R) Xeon(R) Silver 4310 CPU @ 2.10GHz\', 1, \'CpuPhysicCores\', \'24\', 1, \'CpuLogicCores\', \'48\', 1, \'Cpufreqkhz\', \'2100000000.000000\', 1, \'Cpucachelinesizebytes\', \'64\', 1, \'Cpuestlastlevelcachesizemb\', \'30.000000\', 1, \'Cpuestgflops\', \'1344.000000\', 1, \'Memorybwestkbps\', \'170.000000\', 1, \'Numanodes\', \'2\', 1, \'Virtualmemkb\', \'137439000000.000000\', 1, \'Pagesizekb\', \'4096\')', false, true)
					ProfileItem('analysisinfo', 0, 0, 0, 0, 0, 'I(9, 1, \'Analysisstatus\', \'valid\', 1, \'Numsupernodes\', \'359504\', 1, \'Factornnz\', \'1100219968\', 1, \'Factorestflops\', \'4602580000000\', 1, \'Fbsestflops\', \'3642070000\', 1, \'Rootfactestflops\', \'55822200000\', 1, \'Rootfbsestflops\', \'15191072\', 1, \'Analysistimesec\', \'62.644000\', 1, \'Analysismemkb\', \'2073500.000000\')', false, true)
					ProfileItem('factorinfo', 0, 0, 0, 0, 0, 'I(4, 1, \'Fatorizationstatus\', \'valid\', 1, \'Factorizationnumcores\', \'5\', 1, \'Factorizationtimesec\', \'87.402800\', 1, \'Factorizationmentotalkb\', \'13901900.000000\')', false, true)
					ProfileItem('fbsinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Fbstatus\', \'valid\', 1, \'Fbstype\', \'fullsolve\', 1, \'Fbsmt\', \'false\', 1, \'Fbsmrhs\', \'false\', 1, \'Fbsnumcores\', \'5\', 1, \'Fbsnumsolvestotal\', \'1\', 1, \'Fbsnumsolves\', \'1\', 1, \'Fbsavgsolvetime1solvesec\', \'16.670700\', 1, \'Fbscputimesec\', \'16.670700\', 1, \'Fbsmemorytotalkb\', \'10719600.000000\')', false, true)
					ProfileItem('solverprofile', 0, 0, 0, 0, 0, 'I(2, 1, \'Maxmemkb\', \'13901881\', 1, \'Maxdiskkb\', \'0\')', false, true)
				$end 'ProfileGroup'
			$end 'ProfileGroup'
			ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
			$begin 'ProfileGroup'
				MajorVer=2025
				MinorVer=2
				Name='Adaptive Pass 2'
				$begin 'StartInfo'
					I(1, 'Frequency', '5.8GHz')
				$end 'StartInfo'
				$begin 'TotalInfo'
					I(0, ' ')
				$end 'TotalInfo'
				GroupOptions=0
				TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
				ProfileItem('Adaptive Refine', 107, 0, 111, 0, 1725804, 'I(1, 2, \'Tetrahedra\', 1403085, false)', true, true)
				ProfileItem(' ', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
				ProfileItem('Simulation Setup ', 25, 0, 26, 0, 3313568, 'I(2, 2, \'Tetrahedra\', 995225, false, 1, \'Disk\', \'0 Bytes\')', true, true)
				ProfileItem('Matrix Assembly', 179, 0, 309, 0, 10092332, 'I(3, 2, \'Tetrahedra\', 995225, false, 2, \'Lumped ports\', 1, false, 1, \'Disk (output)\', \'422 MB\')', true, true)
				ProfileItem('Matrix Solve', 306, 0, 1104, 0, 26416076, 'I(5, 1, \'Type\', \'MCS\', 2, \'Cores\', 5, false, 2, \'Matrix size\', 5414869, false, 3, \'Matrix bandwidth\', 18.4405, \'%5.1f\', 1, \'Disk\', \'20.7 MB\')', true, true)
				ProfileItem('Field Recovery', 29, 0, 48, 0, 26416076, 'I(2, 2, \'Excitations\', 1, false, 1, \'Disk\', \'44 MB\')', true, true)
				ProfileItem('Data Transfer', 24, 0, 5, 0, 250452, 'I(1, 0, \'Adaptive Pass 2\')', true, true)
				$begin 'ProfileGroup'
					MajorVer=2025
					MinorVer=2
					Name='APIPms'
					$begin 'StartInfo'
						I(1, 'Timesinceepock', '1794416149')
					$end 'StartInfo'
					$begin 'TotalInfo'
						I(0, ' ')
					$end 'TotalInfo'
					GroupOptions=16
					TaskDataOptions(Memory=8)
					ProfileItem('solverinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Solvertype\', \'shared_memory\', 1, \'Precision\', \'mixed\', 1, \'Solversymmetry\', \'complex_sym\', 1, \'Matrixdim\', \'5414869\', 1, \'Matrixbw\', \'18.443300\', 1, \'Matrixnnz\', \'99867927\', 1, \'Rootdim\', \'8548\', 1, \'Mathtype\', \'mkl\', 1, \'Mpitasks\', \'1\', 1, \'Threadspertasks\', \'0\')', false, true)
					ProfileItem('sysinfo', 0, 0, 0, 0, 0, 'I(12, 1, \'Os\', \'win\', 1, \'Cpuid\', \'Intel(R) Xeon(R) Silver 4310 CPU @ 2.10GHz\', 1, \'CpuPhysicCores\', \'24\', 1, \'CpuLogicCores\', \'48\', 1, \'Cpufreqkhz\', \'2100000000.000000\', 1, \'Cpucachelinesizebytes\', \'64\', 1, \'Cpuestlastlevelcachesizemb\', \'30.000000\', 1, \'Cpuestgflops\', \'1344.000000\', 1, \'Memorybwestkbps\', \'170.000000\', 1, \'Numanodes\', \'2\', 1, \'Virtualmemkb\', \'137439000000.000000\', 1, \'Pagesizekb\', \'4096\')', false, true)
					ProfileItem('analysisinfo', 0, 0, 0, 0, 0, 'I(9, 1, \'Analysisstatus\', \'valid\', 1, \'Numsupernodes\', \'477472\', 1, \'Factornnz\', \'1873211505\', 1, \'Factorestflops\', \'11438800000000\', 1, \'Fbsestflops\', \'6277160000\', 1, \'Rootfactestflops\', \'208196000000\', 1, \'Rootfbsestflops\', \'36534152\', 1, \'Analysistimesec\', \'96.193400\', 1, \'Analysismemkb\', \'2918740.000000\')', false, true)
					ProfileItem('factorinfo', 0, 0, 0, 0, 0, 'I(4, 1, \'Fatorizationstatus\', \'valid\', 1, \'Factorizationnumcores\', \'5\', 1, \'Factorizationtimesec\', \'188.055000\', 1, \'Factorizationmentotalkb\', \'26028100.000000\')', false, true)
					ProfileItem('fbsinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Fbstatus\', \'valid\', 1, \'Fbstype\', \'fullsolve\', 1, \'Fbsmt\', \'false\', 1, \'Fbsmrhs\', \'false\', 1, \'Fbsnumcores\', \'5\', 1, \'Fbsnumsolvestotal\', \'1\', 1, \'Fbsnumsolves\', \'1\', 1, \'Fbsavgsolvetime1solvesec\', \'13.923000\', 1, \'Fbscputimesec\', \'13.923000\', 1, \'Fbsmemorytotalkb\', \'18929200.000000\')', false, true)
					ProfileItem('solverprofile', 0, 0, 0, 0, 0, 'I(2, 1, \'Maxmemkb\', \'26028136\', 1, \'Maxdiskkb\', \'0\')', false, true)
				$end 'ProfileGroup'
				ProfileFootnote('I(1, 3, \'Max Mag. Delta S\', 0.227302, \'%.5f\')', 0)
			$end 'ProfileGroup'
		$end 'ProfileGroup'
	$end 'ProfileGroup'
$end 'Profile'
