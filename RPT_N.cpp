
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<math.h>
#include <cstdio>
#define _CRT_SECURE_NO_WARNINGS
#define pcm 10        //对于蒙卡程序，有自带的统计偏差，因此很难做到两个模型的k完全相等，
//所以设置3pcm作为偏差限，两个模型k偏差在3pcm以内，认为等效
#define inf 100000
#define cyc_max 20   //循环上限
#define NUM 3672    //颗粒数
#define pi 3.1415926
#define h 1
#define avo 6.02E+23  //阿伏伽德罗常数
#define barn 1.0E+24
#define FULL_RATE 0.40//设置填充率

//FILE* fp1;
FILE* fp2;
//char *RPT_filein;
//char *RPT_fileout;
double RPT_R;
double R1, R2, R3;
double deltK1, deltK2, deltK3;
double K_fcm;
double bias_fcm;
double K_rpt1, K_rpt2, K_rpt3;
double bias_rpt1, bias_rpt2, bias_rpt3;
double ZL_U5;
double ZL_U8;
double ZL_N14;
double ZL_Si;
double ZL_Si28;
double ZL_Si29;
double ZL_Si30;
double ZL_C;
double r_fuel;
double r_buf;
double r_pyc1;
double r_sic;
double r_pyc2;
double Rou_h2o;
double nd_u8;
double nd_u5;
double nd_n14;
double nd_buf;
double nd_pyc;
double nd_sic;
double nd_h;
double nd_o;
double L;
double r_fcm;
double r_gap;
double r_clad;
char run_command[100];
char none[500];

//读取输入卡信息
void Read_para()
{
	FILE* fp;
	fp = fopen("Parameters", "r");

	none[100] = fscanf(fp, "%s", none);

	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%lf", &r_fuel);//fuel    0.0380
	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%lf", &r_buf);//buffer  0.0400
	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%lf", &r_pyc1);//PyC1     0.0435
	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%lf", &r_sic);//SiC     0.0470
	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%lf", &r_pyc2);//PyC2     0.0490

	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%lf", &r_fcm);//fuel  0.600
	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%lf", &r_gap);//gap  0.606
	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%lf", &r_clad);//clad 0.691
	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%lf", &L);//water  0.824
	printf("L=%f\n", L);

	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%lf", &nd_u5);//92235.03c   2.2519E-03
	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%lf", &nd_u8);//92238.03c   3.1984E-02
	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%lf", &nd_n14);// 7014.03c   3.4236E-02
	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%lf", &nd_buf);//buffer  5.2675E-2
	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%lf", &nd_pyc);//pyc   9.5317E-2
	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%lf", &nd_sic);//SiC_SiC   4.7859e-2
	printf("nd_sic=%E\n", nd_sic);
	none[100] = fscanf(fp, "%s", none);
	none[100] = fscanf(fp, "%lf", &Rou_h2o);//Rou_h2o  1.0
	fclose(fp);
	fp = NULL;
	nd_o = Rou_h2o * avo / 18.0 / barn;
	nd_h = 2 * nd_o;

	printf("---------------数据读取工作完成---------------\n");

}

//均匀化换算
void CalZL()
{
	/*r_fuel = 0.035;
	r_buf = 0.040;
	r_pyc1 = 0.0435;
	r_sic = 0.047;
	r_pyc2 = 0.049;*/


	//ZL_U5 = NUM * 4.0 * pi * r_fuel * r_fuel * r_fuel * nd_u5 / 3.0;//计算U235的核子数
	//ZL_U8 = NUM * 4.0 * pi * r_fuel * r_fuel * r_fuel * nd_u8 / 3.0;//计算U238的核子数
	//ZL_N14 = NUM * 4.0 * pi * r_fuel * r_fuel * r_fuel * nd_n14 / 3.0;//计算N14的核子数
	//ZL_Si = (NUM * 4.0 * pi * (r_sic * r_sic * r_sic - r_pyc1 * r_pyc1 * r_pyc1) / 3.0 + (pi * r_fcm * r_fcm * h - NUM * 4.0 * pi * r_pyc2 * r_pyc2 * r_pyc2 / 3.0)) * nd_sic;
	/*计算SiC的原子数由两部分组成
			颗粒内的SiC层             										基质里面抛开燃料颗粒的SiC			*/
			//ZL_C = ZL_Si + (NUM * 4.0 * pi * (r_pyc1 * r_pyc1 * r_pyc1 - r_buf * r_buf * r_buf) / 3.0 + NUM * 4.0 * pi * (r_pyc2 * r_pyc2 * r_pyc2 - r_sic * r_sic * r_sic) / 3.0) * nd_pyc;
			//ZL_C = ZL_C + (NUM * 4.0 * pi * (r_buf * r_buf * r_buf - r_fuel * r_fuel * r_fuel) / 3.0) * nd_buf;//计算所有的C

	ZL_U5 = FULL_RATE * pi * r_fcm * r_fcm * r_fuel * r_fuel * r_fuel * nd_u5 / (r_pyc2 * r_pyc2 * r_pyc2);
	ZL_U8 = FULL_RATE * pi * r_fcm * r_fcm * r_fuel * r_fuel * r_fuel * nd_u8 / (r_pyc2 * r_pyc2 * r_pyc2);
	ZL_N14 = FULL_RATE * pi * r_fcm * r_fcm * r_fuel * r_fuel * r_fuel * nd_n14 / (r_pyc2 * r_pyc2 * r_pyc2);
	ZL_Si = FULL_RATE * pi * r_fcm * r_fcm * (r_sic * r_sic * r_sic - r_pyc1 * r_pyc1 * r_pyc1) / (r_pyc2 * r_pyc2 * r_pyc2) * nd_sic + (pi * r_fcm * r_fcm) * (1 - FULL_RATE) * nd_sic;
	ZL_C = ZL_Si + FULL_RATE * pi * r_fcm * r_fcm * (r_pyc1 * r_pyc1 * r_pyc1 - r_buf * r_buf * r_buf + r_pyc2 * r_pyc2 * r_pyc2 - r_sic * r_sic * r_sic) / (r_pyc2 * r_pyc2 * r_pyc2) * nd_pyc;
	ZL_C = ZL_C + FULL_RATE * pi * r_fcm * r_fcm * (r_buf * r_buf * r_buf - r_fuel * r_fuel * r_fuel) / (r_pyc2 * r_pyc2 * r_pyc2) * nd_buf;
	printf("ZL_Si=%E\n", ZL_Si);
	printf("ZL_C=%E\n", ZL_C);
	printf("---------------均匀化工作完成---------------\n");
}
//颗粒模型输入卡自动生成

void Alpha_CARD()
{
	FILE* fp;
	fp = fopen("Input.inp", "w");
	fprintf(fp, "\n");
	fprintf(fp, "!CASEID\n");
	fprintf(fp, "CASEID_BLOCK\n");
	fprintf(fp, "CASEID  quartCore2g\n");
	fprintf(fp, "If_Resonance 1\n");
	fprintf(fp, "ProblemType  0		!0--eigenvalue problem, 1--fixed source problem\n");
	fprintf(fp, "\n");
	fprintf(fp, "MATERIAL_BLOCK\n");
	fprintf(fp, "DH_FCMZ_ANE1.inp\n");
	fprintf(fp, "\n");
	fprintf(fp, "!Geometry Information\n");
	fprintf(fp, "\n");
	fprintf(fp, "GEOMETRY_BLOCK\n");
	fprintf(fp, "1_Parallel_Cylinders\n");
	fprintf(fp, "!1 - the information of cylindrical surfaces(PARA_CYLINDER)\n");
	fprintf(fp, "Number_of_Cylinders 3\n");
	fprintf(fp, "!id   surfaceType   boundaryType   center_x   center_y   center_z   radius\n");
	fprintf(fp, "0         6             3             0.        0.          0.      0.6000\n");
	fprintf(fp, "1         6             3             0.        0.          0.      0.6060\n");
	fprintf(fp, "2         6             3             0.        0.          0.      0.6910\n");
	fprintf(fp, "\n");
	fprintf(fp, "2_Primitives\n");
	fprintf(fp, "!2 - different type of Primitives\n");
	fprintf(fp, "Number_of_PrimitiveType 4\n");
	fprintf(fp, "!id   is_cr  mat_id  mesh_numX  mesh_numY  num_rings  num_sectors  num_bound_surf  bound_surface_id  halfspace\n");
	fprintf(fp, "0      0      0        0           0         3           1              1                0            -1\n");
	//我可以在燃料卡里面写两个材料，0代表是燃料颗粒，1代表是燃料，我觉得不太行因为每一个半径所对应的燃料是不一样的，再说吧。
	fprintf(fp, "1      0      1        0           0         1           1              2                0             1\n");
	fprintf(fp, "                                                                                         1            -1\n");
	fprintf(fp, "2      0      2        0           0         1           1              2                1             1\n");
	fprintf(fp, "                                                                                         2            -1\n");
	fprintf(fp, "3      0      3        0           0         1           1              1                2             1\n");
	fprintf(fp, "\n");
	fprintf(fp, "3_Cells\n");
	fprintf(fp, "!3 - different type of cells\n");
	fprintf(fp, "Number_of_CellType  1\n");
	fprintf(fp, "!type_id  ray_tracing_module_id width_x  width_y  num_cells  cells_id\n");
	fprintf(fp, "0              0             1.65      1.65       4       0 1 2 3	!UO2\n");
	fprintf(fp, "\n");
	fprintf(fp, "4_Lattices\n");
	fprintf(fp, "!4 - defferent type of lattice\n");
	fprintf(fp, "Number_of_Lattice  1\n");
	fprintf(fp, "!id  num_x  num_y  UO2 ASSEMBLY\n");
	fprintf(fp, "0    1    1\n");
	fprintf(fp, "0\n");
	fprintf(fp, "\n");
	fprintf(fp, "5_Assemblies\n");
	fprintf(fp, "!5 - different type of assembly\n");
	fprintf(fp, "Number_of_AssemblyType 1\n");
	fprintf(fp, "! id  num_z\n");
	fprintf(fp, "  0     1\n");
	fprintf(fp, "                    lattices      0\n");
	fprintf(fp, "                     width_z     1.0\n");
	fprintf(fp, "                    Sublayers     1\n");
	fprintf(fp, "                    ifaverage     1\n");
	fprintf(fp, "                   sub_width_z\n");
	fprintf(fp, "\n");
	fprintf(fp, "6_Core\n");
	fprintf(fp, "!6 - core geometry information\n");
	fprintf(fp, "!num_x  num_y\n");
	fprintf(fp, "1 1\n");
	fprintf(fp, "0\n");
	fprintf(fp, "\n");
	fprintf(fp, "8_Ray_Tracing_Module\n");
	fprintf(fp, "!8 - specify the dimensions of the ray tracing module\n");
	fprintf(fp, "!mod_dim   x     y     z\n");
	fprintf(fp, "         1.65   1.65   0.0\n");
	fprintf(fp, "\n");
	fprintf(fp, "COUPLE_BLOCK\n");
	fprintf(fp, "rated_power\n");
	fprintf(fp, "rated_flow\n");
	fprintf(fp, "tinlet\n");
	fprintf(fp, "core_power\n");
	fprintf(fp, "\n");
	fprintf(fp, "OPTION_BLOCK\n");
	fprintf(fp, "!Boundary Condition\n");
	fprintf(fp, "!0 - VACUUM, 1 - REFLECTIVE, 2 - PERIODIC, 3 - BOUNDARY_NONE, 4 - EXTRAPOLATION, 5 - ZERO_SURFACE_FLUX\n");
	fprintf(fp, "!WestBC  NorthBC  EastBC  SouthBC  TopBC  BottomBC\n");
	fprintf(fp, "BC   1        1       1      1      1         1\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "!Maximum number of iteration\n");
	fprintf(fp, "!IterationLimit <NOuterMax> <NGSweep> <NInner>\n");
	fprintf(fp, "MOC_ITER_LIM    1000     1       1		!ITER_LIM Card\n");
	fprintf(fp, "!<MG>  <CG>\n");
	fprintf(fp, "CMFD_ITER_LIM   30    10\n");
	fprintf(fp, "LS_ITER_LIM     30    10\n");
	fprintf(fp, "!ConvergenceCriteria <ECritK> <ECritPhi>\n");
	fprintf(fp, "MOC_CONV_CRIT      1.e-5  1.e-5		!CONV_CRIT Card\n");
	fprintf(fp, "!<MG - k> <MG - F> <CG - k> <CG - F>\n");
	fprintf(fp, "CMFD_CONV_CRIT     5.e-7  5.e-6  2.5E-7 2.5E-6\n");
	fprintf(fp, "!<MG>  <CG>\n");
	fprintf(fp, "LS_CONV_CRIT       1.e-6 1.e-6\n");
	fprintf(fp, "SOR_FACTOR         1.\n");
	fprintf(fp, "!GMRES Parameter\n");
	fprintf(fp, "!<MG>  <CG>\n");
	fprintf(fp, "KrylovDim   1     1\n");
	fprintf(fp, "Precond     1     1\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "!solver <solverType>\n");
	fprintf(fp, "SOLVER      	1	!0 - gauss seidel method, 1 - jacobian method\n");
	fprintf(fp, "!source item\n");
	fprintf(fp, "!SOURCE_TYPE  0		!0 - flat source approximation; 1 - linear source approximation\n");
	fprintf(fp, "!Ray Card\n");
	fprintf(fp, "NumPolarAnglePartition  1\n");
	fprintf(fp, "NumGroupBoundPerPolarAngPart  1\n");
	fprintf(fp, "!< PolQuadType> 0 - TYPolarQuad, 1 - LeonardPolarQuad, 2 - GLPolarQuad, 3 - EqualWeightPolarQuad, 4 - EqualAnglePolarQuad\n");
	fprintf(fp, "!RayTracing <Spacing> <NumAzimuAng> <PolQuadType> <OrderTheta> <GroupBound> <PolarBound> <FSRVolCorrType>\n");
	fprintf(fp, "RAY_PARA        0.01          48           0            3          47          3             0\n");
	fprintf(fp, "!whether or not to use a linear exponential table when evaluating the exponentials in the moc equations\n");
	fprintf(fp, "!linear_exp\n");
	fprintf(fp, "EXP      1		!0 - intrinsic exp(x) function, 1 - linear exponential table\n");
	fprintf(fp, "!K - reference\n");
	fprintf(fp, "KType  1		!0--k inf, 1--k eff\n");
	fprintf(fp, "Ref_K  1.48367\n");
	fprintf(fp, "\n");
	fprintf(fp, "!CMFD\n");
	fprintf(fp, "CMFD   0		!1--on, 0--off\n");
	fprintf(fp, "CMFD_NumGroups     47		!number of groups for group condensation\n");
	fprintf(fp, "!the group uper limit need to be condensited relative to original group structure\n");
	fprintf(fp, "CMFD_GroupBound     0  1  2  3  4  5  6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47\n");
	fprintf(fp, "!info for two - level cmfd\n");
	fprintf(fp, "CMFD_NumCoarseGroups  2\n");
	fprintf(fp, "CMFD_CoarseGroupBound\n");
	fprintf(fp, "0 26 47\n");
	fprintf(fp, "\n");
	fprintf(fp, "!Parallel\n");
	fprintf(fp, "!Number of Processes along x, yand z direction\n");
	fprintf(fp, "!Num_subDom_X  Num_subDom_Y  Num_subDom_Z\n");
	fprintf(fp, "!2  2  1\n");
	fprintf(fp, "1  1  1\n");
	fprintf(fp, "!2 2 1\n");
	fprintf(fp, "\n");
	fprintf(fp, "!Processes ID Map - 2D\n");
	fprintf(fp, "!1  2   3  4\n");
	fprintf(fp, "!Hardware ID in Map - 2D, -1 - no set hardware, 0 - cpu, 1 - gpu\n");
	fprintf(fp, "!0 1\n");
	fprintf(fp, "!0 1\n");
	fprintf(fp, "1\n");
	fprintf(fp, "\n");
	fprintf(fp, "!number of cells in each process along x, yand z direction\n");
	fprintf(fp, "!8  9  9  8  1\n");
	fprintf(fp, "1 1 1\n");
	fprintf(fp, "!8 9 17 1\n");
	fprintf(fp, "\n");
	fprintf(fp, "OUTPUT_BLOCK\n");
	fprintf(fp, "!input information\n");
	fprintf(fp, "Output_inputInfo         0		!output the input information\n");
	fprintf(fp, "Output_XSsLibInfo        0		!output the corss - section libarary's information\n");
	fprintf(fp, "!geometry information\n");
	fprintf(fp, "Output_GeometryInfo      0		!output the cell based ray tracing module information\n");
	fprintf(fp, "!MacroXSs information\n");
	fprintf(fp, "  Output_MacroXSs          1		! output the macroscopic cross sections' information\n");
	fprintf(fp, "!results\n");
	fprintf(fp, "Output_eigenvalue        0		!output the eigenvalue of the problem\n");
	fprintf(fp, "Output_normPinPower      0		!output the normalized pin power\n");
	fprintf(fp, "Output_flux              1		!output the relative neutron flux of FSRs\n");
	fprintf(fp, "!Memory footprint\n");
	fprintf(fp, "Output_memoryFootprint   0		!output the memory footprint change in stact\n");
	fprintf(fp, "!CMFD related\n");
	fprintf(fp, "OUTPUT_CMFDGroupConstant 0		!output the group consitant generated for CMFD calculation\n");
	fprintf(fp, "!based on the solution of 2D MOC\n");
	fprintf(fp, "!generate visualization file\n");
	fprintf(fp, "Gene_VisualFile          0\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "!number of threads per MPI process\n");
	fprintf(fp, "NumOMPThread             1\n");
	fprintf(fp, "!CPU transport sweep function set\n");
	fprintf(fp, "CPU_TSFunctionType       2\n");
	fprintf(fp, "!0 Jacobi_MOCSolver_2D_expTab_recExp\n");
	fprintf(fp, "!1 Jacobi_MOCSolver_2D_diam_expTable\n");
	fprintf(fp, "!2 Jacobi_MOCSolver_2D_expFunc_recExp\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "!number of stream\n");
	fprintf(fp, "GPUNumStreams            3\n");
	fprintf(fp, "!set the block dimension\n");
	fprintf(fp, "GPUBlockDim             128\n");
	fprintf(fp, "!number of GPUs per node\n");
	fprintf(fp, "NumGPUsPerNode           1\n");
	fprintf(fp, "!set GPU transport solver type\n");
	fprintf(fp, "GPU_TSKernelType        3\n");
	fprintf(fp, "!0 GPUTS_Group_ExpInte_OneDir_AngFReg\n");
	fprintf(fp, "!1 GPUTS_Group_DDExpInte_OneDir_AngFReg\n");
	fprintf(fp, "!2 GPUTS_Group_ExpFunc_OneDir_AngFReg\n");
	fprintf(fp, "!3 GPUTS_Group_DD_OneDir_AngFReg\n");
	fclose(fp);
	fp = NULL;
	printf("---------------颗粒模型几何卡完成---------------\n");

	fp = fopen("DH_FCMZ_ANE1.inp", "w");
	fprintf(fp, "!0: N densities for mat type; 1: N densities for every fsr; 2: mat map; 100: DH1; 101: DH2\n");
	fprintf(fp, "INPUT_TYPE 100\n");
	fprintf(fp, "!0 : Helios; 1: SHEM_WIMS; 2: UFG_MOC\n");
	fprintf(fp, "LIB_TYPE 1\n");
	fprintf(fp, "LIB_PATH C:\\Users\\Public\\Alpha\\Library_MC\\\n");
	fprintf(fp, "!0 : fsr // 1: ring;\n");
	fprintf(fp, "RES_TYPE    0\n");
	fprintf(fp, "!0 : matrix // 1: index;\n");
	fprintf(fp, "SCATT_TYPE    0\n");
	fprintf(fp, "!0 : no iteration // 1: iteration  // 2: SHEM_1G // 3: SHEM_1G_CAT\n");
	fprintf(fp, "INTER_TYPE  2\n");
	fprintf(fp, "\n");
	fprintf(fp, "MAT_NUMBER 8\n");
	fprintf(fp, "\n");
	fprintf(fp, "!UO2 Fuel 10 %\n");
	fprintf(fp, "!Id   TYPE   Tem  NUC_N  NUC_ID     Density\n");
	fprintf(fp, "  0    0     293   145      91231   1.0000E-21\n");
	fprintf(fp, "                            91233   1.0000E-21\n");
	fprintf(fp, "                            92232   1.0000E-21\n");
	fprintf(fp, "                            92233   1.0000E-21\n");
	fprintf(fp, "                            92234   1.0000E-21\n");
	fprintf(fp, "                            92235   %.4E\n", nd_u5);
	fprintf(fp, "                            92236   1.0000E-21\n");
	fprintf(fp, "                            92237   1.0000E-21\n");
	fprintf(fp, "                            92238   %.4E\n", nd_u8);
	fprintf(fp, "                            93237   1.0000E-21\n");
	fprintf(fp, "                            93238   1.0000E-21\n");
	fprintf(fp, "                            93239   1.0000E-21\n");
	fprintf(fp, "                            94238   1.0000E-21\n");
	fprintf(fp, "                            94239   1.0000E-21\n");
	fprintf(fp, "                            94240   1.0000E-21\n");
	fprintf(fp, "                            94241   1.0000E-21\n");
	fprintf(fp, "                            94242   1.0000E-21\n");
	fprintf(fp, "                            95241   1.0000E-21\n");
	fprintf(fp, "                            95342   1.0000E-21\n");
	fprintf(fp, "                            95243   1.0000E-21\n");
	fprintf(fp, "                            96242   1.0000E-21\n");
	fprintf(fp, "                            96243   1.0000E-21\n");
	fprintf(fp, "                            96244   1.0000E-21\n");
	fprintf(fp, "                            96245   1.0000E-21\n");
	fprintf(fp, "                            96246   1.0000E-21\n");
	fprintf(fp, "                            35581   1.0000E-21\n");
	fprintf(fp, "                            36582   1.0000E-21\n");
	fprintf(fp, "                            36583   1.0000E-21\n");
	fprintf(fp, "                            36584   1.0000E-21\n");
	fprintf(fp, "                            36585   1.0000E-21\n");
	fprintf(fp, "                            36586   1.0000E-21\n");
	fprintf(fp, "                            38589   1.0000E-21\n");
	fprintf(fp, "                            38590   1.0000E-21\n");
	fprintf(fp, "                            39589   1.0000E-21\n");
	fprintf(fp, "                            39591   1.0000E-21\n");
	fprintf(fp, "                            40591   1.0000E-21\n");
	fprintf(fp, "                            40593   1.0000E-21\n");
	fprintf(fp, "                            40595   1.0000E-21\n");
	fprintf(fp, "                            40596   1.0000E-21\n");
	fprintf(fp, "                            41595   1.0000E-21\n");
	fprintf(fp, "                            42595   1.0000E-21\n");
	fprintf(fp, "                            42596   1.0000E-21\n");
	fprintf(fp, "                            42597   1.0000E-21\n");
	fprintf(fp, "                            42598   1.0000E-21\n");
	fprintf(fp, "                            42599   1.0000E-21\n");
	fprintf(fp, "                            42600   1.0000E-21\n");
	fprintf(fp, "                            43599   1.0000E-21\n");
	fprintf(fp, "                            44600   1.0000E-21\n");
	fprintf(fp, "                            44601   1.0000E-21\n");
	fprintf(fp, "                            44602   1.0000E-21\n");
	fprintf(fp, "                            44603   1.0000E-21\n");
	fprintf(fp, "                            44604   1.0000E-21\n");
	fprintf(fp, "                            44605   1.0000E-21\n");
	fprintf(fp, "                            44606   1.0000E-21\n");
	fprintf(fp, "                            45603   1.0000E-21\n");
	fprintf(fp, "                            45605   1.0000E-21\n");
	fprintf(fp, "                            46604   1.0000E-21\n");
	fprintf(fp, "                            46605   1.0000E-21\n");
	fprintf(fp, "                            46606   1.0000E-21\n");
	fprintf(fp, "                            46607   1.0000E-21\n");
	fprintf(fp, "                            46608   1.0000E-21\n");
	fprintf(fp, "                            47609   1.0000E-21\n");
	fprintf(fp, "                            47710   1.0000E-21\n");
	fprintf(fp, "                            47611   1.0000E-21\n");
	fprintf(fp, "                            48610   1.0000E-21\n");
	fprintf(fp, "                            48611   1.0000E-21\n");
	fprintf(fp, "                            48613   1.0000E-21\n");
	fprintf(fp, "                            49615   1.0000E-21\n");
	fprintf(fp, "                            51500   1.0000E-21\n");
	fprintf(fp, "                            51625   1.0000E-21\n");
	fprintf(fp, "                            51627   1.0000E-21\n");
	fprintf(fp, "                            52727   1.0000E-21\n");
	fprintf(fp, "                            52729   1.0000E-21\n");
	fprintf(fp, "                            52632   1.0000E-21\n");
	fprintf(fp, "                            53627   1.0000E-21\n");
	fprintf(fp, "                            53629   1.0000E-21\n");
	fprintf(fp, "                            53631   1.0000E-21\n");
	fprintf(fp, "                            53635   1.0000E-21\n");
	fprintf(fp, "                            54628   1.0000E-21\n");
	fprintf(fp, "                            54630   1.0000E-21\n");
	fprintf(fp, "                            54631   1.0000E-21\n");
	fprintf(fp, "                            54632   1.0000E-21\n");
	fprintf(fp, "                            54633   1.0000E-21\n");
	fprintf(fp, "                            54634   1.0000E-21\n");
	fprintf(fp, "                            54635   1.0000E-21\n");
	fprintf(fp, "                            54636   1.0000E-21\n");
	fprintf(fp, "                            55633   1.0000E-21\n");
	fprintf(fp, "                            55634   1.0000E-21\n");
	fprintf(fp, "                            55635   1.0000E-21\n");
	fprintf(fp, "                            55636   1.0000E-21\n");
	fprintf(fp, "                            55637   1.0000E-21\n");
	fprintf(fp, "                            56634   1.0000E-21\n");
	fprintf(fp, "                            56637   1.0000E-21\n");
	fprintf(fp, "                            56640   1.0000E-21\n");
	fprintf(fp, "                            57639   1.0000E-21\n");
	fprintf(fp, "                            57640   1.0000E-21\n");
	fprintf(fp, "                            58640   1.0000E-21\n");
	fprintf(fp, "                            58641   1.0000E-21\n");
	fprintf(fp, "                            58642   1.0000E-21\n");
	fprintf(fp, "                            58643   1.0000E-21\n");
	fprintf(fp, "                            58644   1.0000E-21\n");
	fprintf(fp, "                            59641   1.0000E-21\n");
	fprintf(fp, "                            59643   1.0000E-21\n");
	fprintf(fp, "                            60642   1.0000E-21\n");
	fprintf(fp, "                            60643   1.0000E-21\n");
	fprintf(fp, "                            60644   1.0000E-21\n");
	fprintf(fp, "                            60645   1.0000E-21\n");
	fprintf(fp, "                            60646   1.0000E-21\n");
	fprintf(fp, "                            60647   1.0000E-21\n");
	fprintf(fp, "                            60648   1.0000E-21\n");
	fprintf(fp, "                            60650   1.0000E-21\n");
	fprintf(fp, "                            61647   1.0000E-21\n");
	fprintf(fp, "                            61648   1.0000E-21\n");
	fprintf(fp, "                            61748   1.0000E-21\n");
	fprintf(fp, "                            61649   1.0000E-21\n");
	fprintf(fp, "                            61651   1.0000E-21\n");
	fprintf(fp, "                            62647   1.0000E-21\n");
	fprintf(fp, "                            62648   1.0000E-21\n");
	fprintf(fp, "                            62649   1.0000E-21\n");
	fprintf(fp, "                            62650   1.0000E-21\n");
	fprintf(fp, "                            62651   1.0000E-21\n");
	fprintf(fp, "                            62652   1.0000E-21\n");
	fprintf(fp, "                            62653   1.0000E-21\n");
	fprintf(fp, "                            62654   1.0000E-21\n");
	fprintf(fp, "                            63651   1.0000E-21\n");
	fprintf(fp, "                            63653   1.0000E-21\n");
	fprintf(fp, "                            63654   1.0000E-21\n");
	fprintf(fp, "                            63655   1.0000E-21\n");
	fprintf(fp, "                            63656   1.0000E-21\n");
	fprintf(fp, "                            63657   1.0000E-21\n");
	fprintf(fp, "                            64654   1.0000E-21\n");
	fprintf(fp, "                            64655   1.0000E-21\n");
	fprintf(fp, "                            64656   1.0000E-21\n");
	fprintf(fp, "                            64657   1.0000E-21\n");
	fprintf(fp, "                            64658   1.0000E-21\n");
	fprintf(fp, "                            64660   1.0000E-21\n");
	fprintf(fp, "                            65659   1.0000E-21\n");
	fprintf(fp, "                            65660   1.0000E-21\n");
	fprintf(fp, "                            65661   1.0000E-21\n");
	fprintf(fp, "                            66661   1.0000E-21\n");
	fprintf(fp, "                            66662   1.0000E-21\n");
	fprintf(fp, "                            66663   1.0000E-21\n");
	fprintf(fp, "                            66664   1.0000E-21\n");
	fprintf(fp, "                            67665   1.0000E-21\n");
	fprintf(fp, "                             8001   %.4E\n", nd_n14);
	fprintf(fp, "!UO2 Fuel 14.3 %\n");
	fprintf(fp, "!Id   TYPE   Tem  NUC_N  NUC_ID     Density\n");
	fprintf(fp, "1    0     293    3      92238   5.39455E-03\n");
	fprintf(fp, "                         92235   2.58587E-02\n");
	fprintf(fp, "                          6000    3.12532E-02\n");
	fprintf(fp, "\n");
	fprintf(fp, "!Carbon buffer\n");
	fprintf(fp, "!Id   TYPE   Tem  NUC_N  NUC_ID     Density\n");
	fprintf(fp, "2     1    293    1      60001      %.6E\n", nd_buf);
	fprintf(fp, "\n");
	fprintf(fp, "!Carbon IPyC / OPyC\n");
	fprintf(fp, "!Id   TYPE   Tem  NUC_N  NUC_ID     Density\n");
	fprintf(fp, "3     1    293    1      60001      %.6E\n", nd_pyc);
	fprintf(fp, "\n");
	fprintf(fp, "!SiC\n");
	fprintf(fp, "!Id   TYPE   Tem  NUC_N  NUC_ID     Density\n");
	fprintf(fp, "4     1    293    4      6000       %.6E\n", nd_sic);
	fprintf(fp, "                         14028      %.6E\n", 0.9223 * nd_sic);
	fprintf(fp, "                         14029      %.6E\n", 0.0467 * nd_sic);
	fprintf(fp, "                         14030      %.6E\n", 0.0310 * nd_sic);
	fprintf(fp, "\n");
	fprintf(fp, "!Clad\n");
	fprintf(fp, "!Id   TYPE   Tem  NUC_N  NUC_ID     Density\n");
	fprintf(fp, "   5     1    293    5     26054    1.896163E-03\n");
	fprintf(fp, "                           26056    2.976709E-02\n");
	fprintf(fp, "                           26057    6.876197E-04\n");
	fprintf(fp, "                           24052    3.235087E-02\n");
	fprintf(fp, "                           13027    3.235087E-02\n");
	fprintf(fp, "\n");
	fprintf(fp, "!Moderator\n");
	fprintf(fp, "!Id   TYPE   Tem  NUC_N  NUC_ID     Density\n");
	fprintf(fp, "6     1    293    2     8016      %.6E\n", nd_o);
	fprintf(fp, "                        1001      %.6E\n", nd_h);
	fprintf(fp, "\n");
	fprintf(fp, "!Gap\n");
	fprintf(fp, "!Id   TYPE   Tem  NUC_N  NUC_ID     Density\n");
	fprintf(fp, "7     1    293    2     2004     2.6900E-05\n");
	fprintf(fp, "                        8016     2.6900E-20\n");
	fprintf(fp, "\n");
	fprintf(fp, "TRISO_NUMBER 1\n");
	fprintf(fp, "\n");
	fprintf(fp, "!tri_id   lay_n\n");
	fprintf(fp, "0        5\n");
	fprintf(fp, "  0.0380  0.0400  0.04350  0.0470  0.0490\n");
	fprintf(fp, "   0       2       3       4       3\n");
	fprintf(fp, "   1       1       1       1       1\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "COM_MAT_NUMBER 4\n");
	fprintf(fp, "\n");
	fprintf(fp, "!Fuel\n");
	fprintf(fp, "!Com1_id  matrix   homo ? n_triso   tri_id    Vol_ratio\n");
	fprintf(fp, "    0     4        0        1         0         %f\n", FULL_RATE);
	fprintf(fp, "!Gap\n");
	fprintf(fp, "!Com1_id  matrix   homo ? n_triso   tri_id    Vol_ratio\n");
	fprintf(fp, "       1     0        0        0       7\n");
	fprintf(fp, "\n");
	fprintf(fp, "!Clad\n");
	fprintf(fp, "!Com1_id  matrix   homo ? n_triso   mat\n");
	fprintf(fp, "    2     0        0        0       5\n");
	fprintf(fp, "\n");
	fprintf(fp, "!Moderator\n");
	fprintf(fp, "!Com1_id  matrix   homo ? n_triso   mat\n");
	fprintf(fp, "3     0        0        0       6\n");
	fprintf(fp, "\n");
	//fprintf(fp, "!OUTP_OPTION");
	fprintf(fp, "IF_OUTP_ISOTOPE_DENSITY          0\n");
	fprintf(fp, "IF_OUTP_FIXED_SRC_XS             1\n");
	fprintf(fp, "IF_OUTP_FIXED_SRC_FLUX           1\n");
	fprintf(fp, "IF_OUTP_MAT_XS                   1\n");
	fprintf(fp, "IF_OUTP_RES_INFO                 1\n");
	fprintf(fp, "IF_OUTP_SCATT_INDEX              1\n");
	fprintf(fp, "IF_OUTP_FIX_XS_CHECK             1\n");
	fprintf(fp, "IF_OUTP_FIX_FLUX_CHECK           1\n");
	fprintf(fp, "IF_OUTP_RES_XS                   1\n");
	fclose(fp);
	fp = NULL;
	printf("---------------颗粒模型材料卡工作完成---------------\n");
}

//等效模型输入卡自动生成
void RPT_CARD(double R)
{
	double ND_fuel;
	double ND_U5;
	double ND_U8;
	double ND_N14;
	double ND_Si;
	double ND_C;
	char RPT_filein[500];
	//RPT_filein = (char*)malloc(500);

	ND_U5 = ZL_U5 / (pi * R * R * h);//U235的核子数密度
	ND_U8 = ZL_U8 / (pi * R * R * h);//U238的核子数密度
	ND_N14 = ZL_N14 / (pi * R * R * h);//N14的核子数密度

	ND_Si = (ZL_Si - (r_fcm * r_fcm - R * R) * pi * h * nd_sic) / (pi * R * R * h);//等效半径内的Si的核子数密度
	ND_C = (ZL_C - (r_fcm * r_fcm - R * R) * pi * h * nd_sic) / (pi * R * R * h);//等效半径内的C的核子数密度
	ND_fuel = ND_U5 + ND_U8 + ND_N14 + ND_Si + ND_C;//fuel的核子数密度



	FILE* fp;
	fp = fopen("Input.inp", "w");
	fprintf(fp, "\n");
	fprintf(fp, "!CASEID\n");
	fprintf(fp, "CASEID_BLOCK\n");
	fprintf(fp, "CASEID  quartCore2g\n");
	fprintf(fp, "If_Resonance 1\n");
	fprintf(fp, "ProblemType  0		!0--eigenvalue problem, 1--fixed source problem\n");
	fprintf(fp, "\n");
	fprintf(fp, "MATERIAL_BLOCK\n");
	fprintf(fp, "DH_FCMZ_ANE1.inp\n");
	fprintf(fp, "\n");
	fprintf(fp, "!Geometry Information\n");
	fprintf(fp, "\n");
	fprintf(fp, "GEOMETRY_BLOCK\n");
	fprintf(fp, "1_Parallel_Cylinders\n");
	fprintf(fp, "!1 - the information of cylindrical surfaces(PARA_CYLINDER)\n");
	if (R != r_fcm)

	{
		fprintf(fp, "Number_of_Cylinders 4\n");
		fprintf(fp, "!id   surfaceType   boundaryType   center_x   center_y   center_z   radius\n");
		fprintf(fp, "0         6             3             0.        0.          0.      %f\n", R);
		fprintf(fp, "1         6             3             0.        0.          0.      %.4f\n", r_fcm);
		fprintf(fp, "2         6             3             0.        0.          0.      %.4f\n", r_gap);
		fprintf(fp, "3         6             3             0.        0.          0.      %.4f\n", r_clad);
		fprintf(fp, "\n");
		fprintf(fp, "2_Primitives\n");
		fprintf(fp, "!2 - different type of Primitives\n");
		fprintf(fp, "Number_of_PrimitiveType 5\n");
		//我可以在燃料卡里面写两个材料，0代表是燃料颗粒，1代表是燃料，我觉得不太行因为每一个半径所对应的燃料是不一样的，再说吧
		fprintf(fp, "!id   is_cr  mat_id  mesh_numX  mesh_numY  num_rings  num_sectors  num_bound_surf  bound_surface_id  halfspace\n");
		fprintf(fp, "0      0      0        0           0         3           1              1               0            -1\n");
		fprintf(fp, "1      0      4        0           0         1           1              2               0             1\n");
		fprintf(fp, "                                                                                        1            -1\n");
		fprintf(fp, "2      0      1        0           0         1           1              2               1             1\n");
		fprintf(fp, "                                                                                        2            -1\n");
		fprintf(fp, "3      0      2        0           0         1           1              2               2             1\n");
		fprintf(fp, "                                                                                        3            -1\n");
		fprintf(fp, "4      0      3        0           0         1           1              1               3             1\n");
		fprintf(fp, "\n");
		fprintf(fp, "3_Cells\n");
		fprintf(fp, "!3 - different type of cells\n");
		fprintf(fp, "Number_of_CellType  1\n");
		fprintf(fp, "!type_id  ray_tracing_module_id width_x  width_y  num_cells  cells_id\n");
		fprintf(fp, "0              0             1.65      1.65       5       0 1 2 3 4	!UO2\n");
	}
	else
	{
		fprintf(fp, "Number_of_Cylinders 3\n");
		fprintf(fp, "!id   surfaceType   boundaryType   center_x   center_y   center_z   radius\n");
		fprintf(fp, "0         6             3             0.        0.          0.      %f\n", R);
		fprintf(fp, "1         6             3             0.        0.          0.      %.4f\n", r_gap);
		fprintf(fp, "2         6             3             0.        0.          0.      %.4f\n", r_clad);
		fprintf(fp, "\n");
		fprintf(fp, "2_Primitives\n");
		fprintf(fp, "!2 - different type of Primitives\n");
		fprintf(fp, "Number_of_PrimitiveType 4\n");
		//我可以在燃料卡里面写两个材料，0代表是燃料颗粒，1代表是燃料，我觉得不太行因为每一个半径所对应的燃料是不一样的，再说吧
		fprintf(fp, "!id   is_cr  mat_id  mesh_numX  mesh_numY  num_rings  num_sectors  num_bound_surf  bound_surface_id  halfspace\n");
		fprintf(fp, "0      0      0        0           0         3           1              1               0            -1\n");
		fprintf(fp, "1      0      1        0           0         1           1              2               0             1\n");
		fprintf(fp, "                                                                                        1            -1\n");
		fprintf(fp, "2      0      2        0           0         1           1              2               1             1\n");
		fprintf(fp, "                                                                                        2            -1\n");
		fprintf(fp, "3      0      3        0           0         1           1              1               2             1\n");
		fprintf(fp, "\n");
		fprintf(fp, "  3_Cells\n");
		fprintf(fp, "! 3 - different type of cells\n");
		fprintf(fp, "  Number_of_CellType  1\n");
		fprintf(fp, "! type_id  ray_tracing_module_id width_x  width_y  num_cells  cells_id\n");
		fprintf(fp, "     0              0             1.65      1.65       4       0 1 2 3	!UO2\n");
	}
	fprintf(fp, "\n");
	fprintf(fp, "4_Lattices\n");
	fprintf(fp, "!4 - defferent type of lattice\n");
	fprintf(fp, "Number_of_Lattice  1\n");
	fprintf(fp, "!id  num_x  num_y  UO2 ASSEMBLY\n");
	fprintf(fp, "0    1    1\n");
	fprintf(fp, "0\n");
	fprintf(fp, "\n");
	fprintf(fp, "5_Assemblies\n");
	fprintf(fp, "!5 - different type of assembly\n");
	fprintf(fp, "Number_of_AssemblyType 1\n");
	fprintf(fp, "!id  num_z\n");
	fprintf(fp, "   0     1\n");
	fprintf(fp, "              lattices      0\n");
	fprintf(fp, "              width_z     1.0\n");
	fprintf(fp, "              Sublayers     1\n");
	fprintf(fp, "              ifaverage     1\n");
	fprintf(fp, "              sub_width_z\n");
	fprintf(fp, "\n");
	fprintf(fp, "6_Core\n");
	fprintf(fp, "!6 - core geometry information\n");
	fprintf(fp, "!num_x  num_y\n");
	fprintf(fp, "1 1\n");
	fprintf(fp, "0\n");
	fprintf(fp, "\n");
	fprintf(fp, "8_Ray_Tracing_Module\n");
	fprintf(fp, "!8 - specify the dimensions of the ray tracing module\n");
	fprintf(fp, "!mod_dim   x     y     z\n");
	fprintf(fp, "         1.65   1.65   0.0\n");
	fprintf(fp, "\n");
	fprintf(fp, "COUPLE_BLOCK\n");
	fprintf(fp, "rated_power\n");
	fprintf(fp, "rated_flow\n");
	fprintf(fp, "tinlet\n");
	fprintf(fp, "core_power\n");
	fprintf(fp, "\n");
	fprintf(fp, "OPTION_BLOCK\n");
	fprintf(fp, "!Boundary Condition\n");
	fprintf(fp, "!0 - VACUUM, 1 - REFLECTIVE, 2 - PERIODIC, 3 - BOUNDARY_NONE, 4 - EXTRAPOLATION, 5 - ZERO_SURFACE_FLUX\n");
	fprintf(fp, "!WestBC  NorthBC  EastBC  SouthBC  TopBC  BottomBC\n");
	fprintf(fp, "BC   1        1       1      1      1         1\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "!Maximum number of iteration\n");
	fprintf(fp, "!IterationLimit <NOuterMax> <NGSweep> <NInner>\n");
	fprintf(fp, "MOC_ITER_LIM    1000     1       1		!ITER_LIM Card\n");
	fprintf(fp, "!<MG>  <CG>\n");
	fprintf(fp, "CMFD_ITER_LIM   30    10\n");
	fprintf(fp, "LS_ITER_LIM     30    10\n");
	fprintf(fp, "!ConvergenceCriteria <ECritK> <ECritPhi>\n");
	fprintf(fp, "MOC_CONV_CRIT      1.e-5  1.e-5		!CONV_CRIT Card\n");
	fprintf(fp, "!<MG - k> <MG - F> <CG - k> <CG - F>\n");
	fprintf(fp, "CMFD_CONV_CRIT     5.e-7  5.e-6  2.5E-7 2.5E-6\n");
	fprintf(fp, "!<MG>  <CG>\n");
	fprintf(fp, "LS_CONV_CRIT       1.e-6 1.e-6\n");
	fprintf(fp, "SOR_FACTOR         1.\n");
	fprintf(fp, "!GMRES Parameter\n");
	fprintf(fp, "!          <MG>  <CG>\n");
	fprintf(fp, "KrylovDim   1     1\n");
	fprintf(fp, "Precond     1     1\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "!solver <solverType>\n");
	fprintf(fp, "SOLVER      	1	! 0-gauss seidel method, 1-jacobian method\n");
	fprintf(fp, "!source item\n");
	fprintf(fp, "!SOURCE_TYPE   0	! 0-flat source approximation; 1-linear source approximation\n");
	fprintf(fp, "!Ray Card\n");
	fprintf(fp, "NumPolarAnglePartition  1\n");
	fprintf(fp, "NumGroupBoundPerPolarAngPart  1\n");
	fprintf(fp, "!< PolQuadType> 0 - TYPolarQuad, 1 - LeonardPolarQuad, 2 - GLPolarQuad, 3 - EqualWeightPolarQuad, 4 - EqualAnglePolarQuad\n");
	fprintf(fp, "!RayTracing <Spacing> <NumAzimuAng> <PolQuadType> <OrderTheta> <GroupBound> <PolarBound> <FSRVolCorrType>\n");
	fprintf(fp, "RAY_PARA        0.01          64           0            3          47          3             0\n");
	fprintf(fp, "!whether or not to use a linear exponential table when evaluating the exponentials in the moc equations\n");
	fprintf(fp, "!linear_exp\n");
	fprintf(fp, "EXP      1		!0 - intrinsic exp(x) function, 1 - linear exponential table\n");
	fprintf(fp, "!K - reference\n");
	fprintf(fp, "KType  1		!0--k inf, 1--k eff\n");
	fprintf(fp, "Ref_K  1.48367\n");
	fprintf(fp, "\n");
	fprintf(fp, "!CMFD\n");
	fprintf(fp, "CMFD   0		!1--on, 0--off\n");
	fprintf(fp, "CMFD_NumGroups     47		!number of groups for group condensation\n");
	fprintf(fp, "!the group uper limit need to be condensited relative to original group structure\n");
	fprintf(fp, "CMFD_GroupBound     0  1  2  3  4  5  6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47\n");
	fprintf(fp, "!info for two - level cmfd\n");
	fprintf(fp, "CMFD_NumCoarseGroups  2\n");
	fprintf(fp, "CMFD_CoarseGroupBound\n");
	fprintf(fp, "0 26 47\n");
	fprintf(fp, "\n");
	fprintf(fp, "!Parallel\n");
	fprintf(fp, "!Number of Processes along x, yand z direction\n");
	fprintf(fp, "!Num_subDom_X  Num_subDom_Y  Num_subDom_Z\n");
	fprintf(fp, "!2  2  1\n");
	fprintf(fp, "1  1  1\n");
	fprintf(fp, "!2 2 1\n");
	fprintf(fp, "\n");
	fprintf(fp, "!Processes ID Map - 2D\n");
	fprintf(fp, "!1  2   3  4\n");
	fprintf(fp, "!Hardware ID in Map - 2D, -1 - no set hardware, 0 - cpu, 1 - gpu\n");
	fprintf(fp, "!0 1\n");
	fprintf(fp, "!0 1\n");
	fprintf(fp, "1\n");
	fprintf(fp, "\n");
	fprintf(fp, "!number of cells in each process along x, yand z direction\n");
	fprintf(fp, "!8  9  9  8  1\n");
	fprintf(fp, "1 1 1\n");
	fprintf(fp, "!8 9 17 1\n");
	fprintf(fp, "\n");
	fprintf(fp, "OUTPUT_BLOCK\n");
	fprintf(fp, "!input information\n");
	fprintf(fp, "Output_inputInfo         0		!output the input information\n");
	fprintf(fp, "Output_XSsLibInfo        0		!output the corss - section libarary's information\n");
	fprintf(fp, "!geometry information\n");
	fprintf(fp, "Output_GeometryInfo      0		!output the cell based ray tracing module information\n");
	fprintf(fp, "!MacroXSs information\n");
	fprintf(fp, "Output_MacroXSs          0		!output the macroscopic cross sections' information\n");
	fprintf(fp, "!results\n");
	fprintf(fp, "Output_eigenvalue        0		!output the eigenvalue of the problem\n");
	fprintf(fp, "Output_normPinPower      0		!output the normalized pin power\n");
	fprintf(fp, "Output_flux              1		!output the relative neutron flux of FSRs\n");
	fprintf(fp, "!Memory footprint\n");
	fprintf(fp, "Output_memoryFootprint   0		!output the memory footprint change in stact\n");
	fprintf(fp, "!CMFD related\n");
	fprintf(fp, "OUTPUT_CMFDGroupConstant 0		!output the group consitant generated for CMFD calculation\n");
	fprintf(fp, "!based on the solution of 2D MOC\n");
	fprintf(fp, "!generate visualization file\n");
	fprintf(fp, "Gene_VisualFile          0\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "!number of threads per MPI process\n");
	fprintf(fp, "NumOMPThread             1\n");
	fprintf(fp, "!CPU transport sweep function set\n");
	fprintf(fp, "CPU_TSFunctionType       2\n");
	fprintf(fp, "!0 Jacobi_MOCSolver_2D_expTab_recExp\n");
	fprintf(fp, "!1 Jacobi_MOCSolver_2D_diam_expTable\n");
	fprintf(fp, "!2 Jacobi_MOCSolver_2D_expFunc_recExp\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "!number of stream\n");
	fprintf(fp, "GPUNumStreams            3\n");
	fprintf(fp, "!set the block dimension\n");
	fprintf(fp, "GPUBlockDim             128\n");
	fprintf(fp, "!number of GPUs per node\n");
	fprintf(fp, "NumGPUsPerNode           1\n");
	fprintf(fp, "!set GPU transport solver type\n");
	fprintf(fp, "GPU_TSKernelType        3\n");
	fprintf(fp, "!0 GPUTS_Group_ExpInte_OneDir_AngFReg\n");
	fprintf(fp, "!1 GPUTS_Group_DDExpInte_OneDir_AngFReg\n");
	fprintf(fp, "!2 GPUTS_Group_ExpFunc_OneDir_AngFReg\n");
	fprintf(fp, "!3 GPUTS_Group_DD_OneDir_AngFReg\n");
	fclose(fp);
	fp = NULL;
	printf("---------------RPT几何卡完成---------------\n");



	fp = fopen("DH_FCMZ_ANE1.inp", "w");
	fprintf(fp, "!0: N densities for mat type; 1: N densities for every fsr; 2: mat map\n");
	fprintf(fp, "INPUT_TYPE 0\n");
	fprintf(fp, "!0 : Helios; 1: SHEM_WIMS; 2: UFG_MOC\n");
	fprintf(fp, "LIB_TYPE 1\n");
	fprintf(fp, "LIB_PATH C:\\Users\\Public\\Alpha\\Library_MC\\\n");
	fprintf(fp, "!0: fsr // 1: ring;\n");
	fprintf(fp, "RES_TYPE    0\n");
	fprintf(fp, "!0 : matrix // 1: index;\n");
	fprintf(fp, "SCATT_TYPE    0\n");
	fprintf(fp, "!0 : 408g no iteration // 1: 47g helios iteration  // 2: 47g SHEM_1G // 3: (not use) SHEM_1G_CAT\n");
	fprintf(fp, "INTER_TYPE  2\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "MAT_NUMBER 5\n");
	fprintf(fp, "\n");
	fprintf(fp, "!UO2 Fuel\n");
	fprintf(fp, "!Id   TYPE   Tem  NUC_N  NUC_ID     Density\n");
	fprintf(fp, " 0     0     293    7    92235      %.6E\n", ND_U5);
	fprintf(fp, "                         92238      %.6E\n", ND_U8);
	fprintf(fp, "                         14028      %.6E\n", 0.9223 * ND_Si);
	fprintf(fp, "                         14029      %.6E\n", 0.0467 * ND_Si);
	fprintf(fp, "                         14030      %.6E\n", 0.0310 * ND_Si);
	fprintf(fp, "                         8001       %.6E\n", ND_N14);
	fprintf(fp, "                         6000       %.6E\n", ND_C);
	fprintf(fp, "\n");
	fprintf(fp, "!Gap\n");
	fprintf(fp, "!Id   TYPE   Tem  NUC_N  NUC_ID     Density\n");
	fprintf(fp, " 1     1     293    2     2004      2.69E-5\n");
	fprintf(fp, "                          8016      1.0E-20");
	fprintf(fp, "\n");
	fprintf(fp, "!clad\n");
	fprintf(fp, "!Id   TYPE   Tem  NUC_N  NUC_ID     Density\n");
	fprintf(fp, " 2      0     293   5    26054    1.896163E-03 \n");
	fprintf(fp, "                          26056    2.976709E-02\n");
	fprintf(fp, "                          26057    6.876197E-04\n");
	fprintf(fp, "                          24052    3.235087E-02\n");
	fprintf(fp, "                          13027    3.235087E-02\n");
	fprintf(fp, "\n");
	fprintf(fp, "!Moderator\n");
	fprintf(fp, "!Id   TYPE   Tem  NUC_N  NUC_ID     Density\n");
	fprintf(fp, "3     1      293  2     8016         %.6E\n", nd_o);
	fprintf(fp, "                        1001         %.6E\n", nd_h);
	fprintf(fp, "\n");
	fprintf(fp, "!clad with gap\n");
	fprintf(fp, "!Id   TYPE   Tem  NUC_N  NUC_ID     Density\n");
	fprintf(fp, " 4     1     293    4    6000       %.6E\n", nd_sic);
	fprintf(fp, "                         14028      %.6E\n", 0.9223 * nd_sic);
	fprintf(fp, "                         14029      %.6E\n", 0.0467 * nd_sic);
	fprintf(fp, "                         14030      %.6E\n", 0.0310 * nd_sic);
	fprintf(fp, "\n");
	fprintf(fp, "!OUTP_OPTION\n");
	fprintf(fp, "IF_OUTP_ISOTOPE_DENSITY          1\n");
	fprintf(fp, "IF_OUTP_FIXED_SRC_XS             1\n");
	fprintf(fp, "IF_OUTP_FIXED_SRC_FLUX           1\n");
	fprintf(fp, "IF_OUTP_MAT_XS                   1\n");
	fprintf(fp, "IF_OUTP_RES_INFO                 1\n");
	fprintf(fp, "IF_OUTP_SCATT_INDEX              1\n");
	fprintf(fp, "IF_OUTP_FIX_XS_CHECK             1\n");
	fprintf(fp, "IF_OUTP_FIX_FLUX_CHECK           1\n");
	fprintf(fp, "IF_OUTP_RES_XS                   1\n");
	fprintf(fp, "\n\n\n\n");

	fclose(fp);
	fp = NULL;
	printf("---------------RPT材料卡完成---------------\n");
}

void Copy()
{
	char chere;
	FILE* prf = fopen("K_EFF.out", "r");
	FILE* prw = fopen("KEFF.txt", "w");
	if (NULL == prf)
	{
		perror("open file K_EFF.out");
	}
	if (NULL == prw)
	{
		perror("open file KEFF.txtx");
	}
	while ((chere = fgetc(prf)) != EOF)
	{
		fputc(chere, prw);
	}
	fclose(prf);
	fclose(prw);
	prf = NULL;
	prw = NULL;

}

//读取k
double ReadFCMK(double K)
{

	char none[500];
	char symbol[] = "IMP_KEFF";
	fp2 = fopen("K_EFF.out", "r");
	// Eigenvalue Problem: k= 1.373787    
	none[100] = fscanf(fp2, "%s", none);//Eigenvalue
	none[100] = fscanf(fp2, "%s", none); //Problem
	//none[100] = fscanf(fp2, "%s", none); //:
	//none[100] = fscanf(fp2, "%s", none); //k
	none[100] = fscanf(fp2, "%s", none); //= 
	none[100] = fscanf(fp2, "%lf", &K);//
	fclose(fp2);
	fp2 = NULL;
	return K;
	printf("---------------读取K_EFF完成---------------\n");

}

int main()
{
	FILE* tp;
	tp = fopen("record.m", "a");
	fprintf(tp, "Record the RPT radius and kinf\n\n");
	char run_command[100];
	char RPT_filein[500];
	int i, j, k;
	double Rmin, Rmax;
	Read_para();
	Rmin = 0.6324555 * r_fcm;      //半径下限（这么说就是40%的填充率咯）
	Rmax = r_fcm;                //半径上限
	system("mkdir Alpha");//新建一个文件夹
	system("copy ALPHA-DH_Win64.exe Alpha");
	system("copy burn.inp Alpha");
	system("copy DepthMainLib Alpha");
	system("copy HLS_Lim_Lib Alpha");
	CalZL();
	Alpha_CARD();
	printf("Alpha CARD OK\n");
	system(".\\ALPHA-DH_Win64.exe");  //运行Alpha
	printf("HAPPY GAME!\n");
	K_fcm = ReadFCMK(K_fcm);
	printf("K_fcm=%f\n", K_fcm);
	Copy();
	system("move Input.inp Alpha");
	system("move DH_FCMZ_ANE1.inp Alpha");
	system("move KEFF.txt Alpha");

	R2 = Rmax;
	system("mkdir Rmax");
	system("copy  ALPHA-DH_Win64.exe Rmax");
	system("copy burn.inp Rmax");
	system("copy DepthMainLib Rmax");
	system("copy HLS_Lim_Lib Rmax");
	//system("cd Rmax");
	RPT_CARD(R2);
	//system("pause");
	system(" .\\ALPHA-DH_Win64.exe");  //运行Alpha
	printf("HAPPY GAME!\n");
	K_rpt2 = ReadFCMK(K_rpt2);
	Copy();
	system("move Input.inp Rmax");
	system("move DH_FCMZ_ANE1.inp Rmax");
	system("move KEFF.txt Rmax");
	printf("K2=%f  Kfcm=%f\n", K_rpt2, K_fcm);
	deltK2 = (K_rpt2 - K_fcm) * 100000;
	//system("cd ..");

	R1 = Rmin;
	system("mkdir Rmin");
	system("copy ALPHA-DH_Win64.exe Rmin");
	system("copy burn.inp Rmin");
	system("copy DepthMainLib Rmin");
	system("copy HLS_Lim_Lib Rmin");
	system("cd Rmin");
	RPT_CARD(R1);
	system(" .\\ALPHA-DH_Win64.exe");  //运行Alpha
	printf("HAPPY GAME!\n");
	K_rpt1 = ReadFCMK(K_rpt1);
	Copy();
	system("move Input.inp Rmin");
	system("move DH_FCMZ_ANE1.inp Rmin");
	system("move KEFF.txt Rmin");
	printf("K2=%f  Kfcm=%f\n", K_rpt1, K_fcm);
	deltK1 = (K_rpt1 - K_fcm) * 100000;
	system("cd ..");

	//double R1 = 0.55, R2 = 0.50, R3 = 0.45;
	//double R4 = 0.400;
	////double R[] = {R1,R2,R3,R4,0.59,0.58,0.57,0.56,0.54,0.53,0.52,0.51,0.49,0.48,0.47,0.46,0.44,0.43,0.43,0.41};
	//double R[20] = {};
	//for (k = 0; k < 10; k++)
	//{
	//	R[k] = 0.435 + 0.0001 * k;
	//	R[k + 10] = 0.502 + 0.0001 * k;
	//}
	//for (k = 0; k < 10; k++)
	//{
	//	R[k] = 0.50294 + 0.000001 * k;
	//}
	//j = sizeof(R) / sizeof(R[0]);
	//double Keff[100] = { 0 };

	//printf("我要开始计算我选定的半径了\n");
	//for (i = 0; i < j; i++)
	//{
	//	sprintf(run_command, "mkdir %f", R[i]);
	//	system(run_command);
	//	sprintf(run_command, "copy ALPHA-DH_Win64.exe %f", R[i]);
	//	system(run_command);
	//	sprintf(run_command, "copy burn.inp %f", R[i]);
	//	system(run_command);
	//	sprintf(run_command, "copy DepthMainLib %f", R[i]);
	//	system(run_command);
	//	sprintf(run_command, "copy HLS_Lim_Lib %f", R[i]);
	//	system(run_command);
	//	RPT_CARD(R[i]);
	//	system(" .\\ALPHA-DH_Win64.exe");  //运行Alpha
	//	printf("HAPPY GAME!\n");
	//	Keff[i] = ReadFCMK(Keff[i]);
	//	Copy();
	//	sprintf(run_command, "move Input.inp %f", R[i]);
	//	system(run_command);
	//	sprintf(run_command, "move DH_FCMZ_ANE1.inp %f", R[i]);
	//	system(run_command);
	//	sprintf(run_command, "move KEFF.txt %f", R[i]);
	//	system(run_command);
	//	printf("happy,game!\n");

	//}
	//printf("沙雕，你选的半径已经算完啦！\n");
	//system("pause");

	//搜索半径
	for (i = 0; i < cyc_max; i++)
	{

		printf("K1=%f  K2=%f  Kfcm=%f\n", K_rpt1, K_rpt2, K_fcm);
		R3 = (R2 - R1) * (K_fcm - K_rpt1) / (K_rpt2 - K_rpt1) + R1;

		printf("R3=%f\n", R3);
		RPT_CARD(R3);
		sprintf(RPT_filein, "%f", R3);
		//sprintf(run_command, "sss2 -omp 40 %s", RPT_filein);
		sprintf(run_command, "mkdir %d", i);
		system(run_command);  //创建一个文件夹
		sprintf(run_command, "copy ALPHA-DH_Win64.exe %d", i);
		system(run_command);  //复制alpha执行程序进入文件夹
		sprintf(run_command, "copy burn.inp %d", i);
		system(run_command);//复制燃耗点设置进入文件夹
		sprintf(run_command, "copy DepthMainLib %d", i);
		system(run_command);//复制路径进入文件夹
		sprintf(run_command, "copy HLS_Lim_Lib %d", i);
		system(run_command);//复制路径进入文件夹
		system(" .\\ALPHA-DH_Win64.exe");  //运行Alpha
		K_rpt3 = ReadFCMK(K_rpt3);
		deltK3 = (K_rpt3 - K_fcm) * 100000;
		Copy();
		sprintf(run_command, "move Input.inp %d", i);
		system(run_command);
		sprintf(run_command, "move DH_FCMZ_ANE1.inp %d", i);
		system(run_command);
		sprintf(run_command, "move KEFF.txt %d", i);
		system(run_command);
		printf("Happy,Game!\n");
		printf("deltK3=%.2f\n", deltK3);
		if ((deltK3 >= (-1.0 * pcm)) && (deltK3 <= pcm))
		{
			i = 100;
		}

		if (deltK3 > pcm)
		{
			R1 = R3;
			K_rpt1 = K_rpt3;
		}
		else if (deltK3 < (-1.0 * pcm))
		{
			R2 = R3;
			K_rpt2 = K_rpt3;
		}


	}
	fclose(tp);
	printf("Finish searching RPT Radius!\n");
	printf("The Radius of RPT method : %f\n", R3);
	printf("The kinf of RPT method : %f \n", K_rpt3);
	printf("The kinf of FCM fuel :   %f \n", K_fcm);

}