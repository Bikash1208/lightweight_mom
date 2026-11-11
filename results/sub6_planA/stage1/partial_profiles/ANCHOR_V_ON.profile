$begin 'Profile'
	$begin 'ProfileGroup'
		MajorVer=2025
		MinorVer=2
		Name='Solution Process'
		$begin 'StartInfo'
			I(1, 'Start Time', '11/11/2026 21:47:17')
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
			ProfileItem('Machine', 0, 0, 0, 0, 0, 'I(6, 1, \'Name\', \'AEM-LAB1\', 1, \'Memory\', \'256 GB\', 3, \'RAM Limit\', 90, \'%f%%\', 2, \'Tasks\', 1, false, 2, \'Cores\', 5, false, 1, \'Free Disk Space\', \'311 GB\')', false, true)
		$end 'ProfileGroup'
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 1, \'Allow off core\', \'True\')', false, true)
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 1, \'Solution Basis Order\', \'1\')', false, true)
		ProfileItem('Design Validation', 0, 0, 0, 0, 0, 'I(1, 0, \'Elapsed time : 00:00:00 , HFSS ComEngine Memory : 189 M\')', false, true)
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'Perform full validations with standard port validations\')', false, true)
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
		$begin 'ProfileGroup'
			MajorVer=2025
			MinorVer=2
			Name='Initial Meshing'
			$begin 'StartInfo'
				I(1, 'Time', '11/11/2026 21:47:19')
			$end 'StartInfo'
			$begin 'TotalInfo'
				I(1, 'Elapsed Time', '00:18:16')
			$end 'TotalInfo'
			GroupOptions=4
			TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
			ProfileItem('Mesh', 603, 0, 808, 0, 5466000, 'I(3, 2, \'Tetrahedra\', 4254476, false, 1, \'Type\', \'TAU\', 2, \'Cores\', 5, false)', true, true)
			ProfileItem('Coarsen', 361, 0, 361, 0, 5466000, 'I(1, 2, \'Tetrahedra\', 1178785, false)', true, true)
			ProfileItem('Lambda Refine', 25, 0, 29, 0, 1147428, 'I(1, 2, \'Tetrahedra\', 1200262, false)', true, true)
			ProfileItem('Simulation Setup', 23, 0, 23, 0, 2811492, 'I(2, 2, \'Tetrahedra\', 788278, false, 1, \'Disk\', \'0 Bytes\')', true, true)
			ProfileItem('Port Adapt', 14, 0, 14, 0, 2830208, 'I(2, 2, \'Tetrahedra\', 788278, false, 1, \'Disk\', \'76.1 KB\')', true, true)
			ProfileItem('Port Refine', 58, 0, 61, 0, 1495296, 'I(1, 2, \'Tetrahedra\', 1196263, false)', true, true)
		$end 'ProfileGroup'
		$begin 'ProfileGroup'
			MajorVer=2025
			MinorVer=2
			Name='Adaptive Meshing'
			$begin 'StartInfo'
				I(1, 'Time', '11/11/2026 22:05:36')
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
				ProfileItem('Simulation Setup ', 19, 0, 20, 0, 2804920, 'I(2, 2, \'Tetrahedra\', 789019, false, 1, \'Disk\', \'185 KB\')', true, true)
				ProfileItem('Matrix Assembly', 97, 0, 166, 0, 7712848, 'I(3, 2, \'Tetrahedra\', 789019, false, 2, \'Lumped ports\', 1, false, 1, \'Disk (output)\', \'305 MB\')', true, true)
				ProfileItem('Matrix Solve', 165, 0, 532, 0, 16804112, 'I(5, 1, \'Type\', \'MCS\', 2, \'Cores\', 5, false, 2, \'Matrix size\', 4114809, false, 3, \'Matrix bandwidth\', 17.455, \'%5.1f\', 1, \'Disk\', \'15.7 MB\')', true, true)
				ProfileItem('Field Recovery', 17, 0, 33, 0, 16804112, 'I(2, 2, \'Excitations\', 1, false, 1, \'Disk\', \'99.4 MB\')', true, true)
				ProfileItem('Data Transfer', 2, 0, 2, 0, 210440, 'I(1, 0, \'Adaptive Pass 1\')', true, true)
				$begin 'ProfileGroup'
					MajorVer=2025
					MinorVer=2
					Name='APIPms'
					$begin 'StartInfo'
						I(1, 'Timesinceepock', '1794415220')
					$end 'StartInfo'
					$begin 'TotalInfo'
						I(0, ' ')
					$end 'TotalInfo'
					GroupOptions=16
					TaskDataOptions(Memory=8)
					ProfileItem('solverinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Solvertype\', \'shared_memory\', 1, \'Precision\', \'mixed\', 1, \'Solversymmetry\', \'complex_sym\', 1, \'Matrixdim\', \'4114809\', 1, \'Matrixbw\', \'17.456499\', 1, \'Matrixnnz\', \'71830145\', 1, \'Rootdim\', \'5740\', 1, \'Mathtype\', \'mkl\', 1, \'Mpitasks\', \'1\', 1, \'Threadspertasks\', \'0\')', false, true)
					ProfileItem('sysinfo', 0, 0, 0, 0, 0, 'I(12, 1, \'Os\', \'win\', 1, \'Cpuid\', \'Intel(R) Xeon(R) Silver 4310 CPU @ 2.10GHz\', 1, \'CpuPhysicCores\', \'24\', 1, \'CpuLogicCores\', \'48\', 1, \'Cpufreqkhz\', \'2100000000.000000\', 1, \'Cpucachelinesizebytes\', \'64\', 1, \'Cpuestlastlevelcachesizemb\', \'30.000000\', 1, \'Cpuestgflops\', \'1344.000000\', 1, \'Memorybwestkbps\', \'170.000000\', 1, \'Numanodes\', \'2\', 1, \'Virtualmemkb\', \'137439000000.000000\', 1, \'Pagesizekb\', \'4096\')', false, true)
					ProfileItem('analysisinfo', 0, 0, 0, 0, 0, 'I(9, 1, \'Analysisstatus\', \'valid\', 1, \'Numsupernodes\', \'365883\', 1, \'Factornnz\', \'1103671343\', 1, \'Factorestflops\', \'4516490000000\', 1, \'Fbsestflops\', \'3651260000\', 1, \'Rootfactestflops\', \'63039800000\', 1, \'Rootfbsestflops\', \'16473800\', 1, \'Analysistimesec\', \'61.134300\', 1, \'Analysismemkb\', \'2106830.000000\')', false, true)
					ProfileItem('factorinfo', 0, 0, 0, 0, 0, 'I(4, 1, \'Fatorizationstatus\', \'valid\', 1, \'Factorizationnumcores\', \'5\', 1, \'Factorizationtimesec\', \'90.253500\', 1, \'Factorizationmentotalkb\', \'13835000.000000\')', false, true)
					ProfileItem('fbsinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Fbstatus\', \'valid\', 1, \'Fbstype\', \'fullsolve\', 1, \'Fbsmt\', \'false\', 1, \'Fbsmrhs\', \'false\', 1, \'Fbsnumcores\', \'5\', 1, \'Fbsnumsolvestotal\', \'1\', 1, \'Fbsnumsolves\', \'1\', 1, \'Fbsavgsolvetime1solvesec\', \'8.850560\', 1, \'Fbscputimesec\', \'8.850560\', 1, \'Fbsmemorytotalkb\', \'10960100.000000\')', false, true)
					ProfileItem('solverprofile', 0, 0, 0, 0, 0, 'I(2, 1, \'Maxmemkb\', \'13835032\', 1, \'Maxdiskkb\', \'0\')', false, true)
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
				ProfileItem('Adaptive Refine', 95, 0, 99, 0, 1815240, 'I(1, 2, \'Tetrahedra\', 1432969, false)', true, true)
				ProfileItem(' ', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
				ProfileItem('Simulation Setup ', 24, 0, 25, 0, 3370032, 'I(2, 2, \'Tetrahedra\', 1014922, false, 1, \'Disk\', \'0 Bytes\')', true, true)
				ProfileItem('Matrix Assembly', 183, 0, 298, 0, 10282840, 'I(3, 2, \'Tetrahedra\', 1014922, false, 2, \'Lumped ports\', 1, false, 1, \'Disk (output)\', \'430 MB\')', true, true)
				ProfileItem('Matrix Solve', 327, 0, 1191, 0, 26570832, 'I(5, 1, \'Type\', \'MCS\', 2, \'Cores\', 5, false, 2, \'Matrix size\', 5521011, false, 3, \'Matrix bandwidth\', 18.4351, \'%5.1f\', 1, \'Disk\', \'21.1 MB\')', true, true)
				ProfileItem('Field Recovery', 38, 0, 61, 0, 26570832, 'I(2, 2, \'Excitations\', 1, false, 1, \'Disk\', \'45 MB\')', true, true)
				ProfileItem('Data Transfer', 32, 0, 8, 0, 210560, 'I(1, 0, \'Adaptive Pass 2\')', true, true)
				$begin 'ProfileGroup'
					MajorVer=2025
					MinorVer=2
					Name='APIPms'
					$begin 'StartInfo'
						I(1, 'Timesinceepock', '1794415878')
					$end 'StartInfo'
					$begin 'TotalInfo'
						I(0, ' ')
					$end 'TotalInfo'
					GroupOptions=16
					TaskDataOptions(Memory=8)
					ProfileItem('solverinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Solvertype\', \'shared_memory\', 1, \'Precision\', \'mixed\', 1, \'Solversymmetry\', \'complex_sym\', 1, \'Matrixdim\', \'5521011\', 1, \'Matrixbw\', \'18.437799\', 1, \'Matrixnnz\', \'101795204\', 1, \'Rootdim\', \'8160\', 1, \'Mathtype\', \'mkl\', 1, \'Mpitasks\', \'1\', 1, \'Threadspertasks\', \'0\')', false, true)
					ProfileItem('sysinfo', 0, 0, 0, 0, 0, 'I(12, 1, \'Os\', \'win\', 1, \'Cpuid\', \'Intel(R) Xeon(R) Silver 4310 CPU @ 2.10GHz\', 1, \'CpuPhysicCores\', \'24\', 1, \'CpuLogicCores\', \'48\', 1, \'Cpufreqkhz\', \'2100000000.000000\', 1, \'Cpucachelinesizebytes\', \'64\', 1, \'Cpuestlastlevelcachesizemb\', \'30.000000\', 1, \'Cpuestgflops\', \'1344.000000\', 1, \'Memorybwestkbps\', \'170.000000\', 1, \'Numanodes\', \'2\', 1, \'Virtualmemkb\', \'137439000000.000000\', 1, \'Pagesizekb\', \'4096\')', false, true)
					ProfileItem('analysisinfo', 0, 0, 0, 0, 0, 'I(9, 1, \'Analysisstatus\', \'valid\', 1, \'Numsupernodes\', \'487093\', 1, \'Factornnz\', \'1892683364\', 1, \'Factorestflops\', \'11173700000000\', 1, \'Fbsestflops\', \'6359540000\', 1, \'Rootfactestflops\', \'181113000000\', 1, \'Rootfbsestflops\', \'33292800\', 1, \'Analysistimesec\', \'94.312700\', 1, \'Analysismemkb\', \'2975140.000000\')', false, true)
					ProfileItem('factorinfo', 0, 0, 0, 0, 0, 'I(4, 1, \'Fatorizationstatus\', \'valid\', 1, \'Factorizationnumcores\', \'5\', 1, \'Factorizationtimesec\', \'209.018000\', 1, \'Factorizationmentotalkb\', \'24752400.000000\')', false, true)
					ProfileItem('fbsinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Fbstatus\', \'valid\', 1, \'Fbstype\', \'fullsolve\', 1, \'Fbsmt\', \'false\', 1, \'Fbsmrhs\', \'false\', 1, \'Fbsnumcores\', \'5\', 1, \'Fbsnumsolvestotal\', \'1\', 1, \'Fbsnumsolves\', \'1\', 1, \'Fbsavgsolvetime1solvesec\', \'16.955900\', 1, \'Fbscputimesec\', \'16.955900\', 1, \'Fbsmemorytotalkb\', \'18943800.000000\')', false, true)
					ProfileItem('solverprofile', 0, 0, 0, 0, 0, 'I(2, 1, \'Maxmemkb\', \'24752388\', 1, \'Maxdiskkb\', \'0\')', false, true)
				$end 'ProfileGroup'
				ProfileFootnote('I(1, 3, \'Max Mag. Delta S\', 0.175612, \'%.5f\')', 0)
			$end 'ProfileGroup'
		$end 'ProfileGroup'
	$end 'ProfileGroup'
$end 'Profile'
