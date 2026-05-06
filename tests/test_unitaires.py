import unittest, os, sys, time
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
import calins as cl


TEST_INPUTS_ROOT = os.path.join(os.path.dirname(__file__), "test_inputs")

sensi_correct_path_1 = os.path.join(TEST_INPUTS_ROOT, "sensi_correct_1.sdf")
sensi_correct_path_2 = os.path.join(TEST_INPUTS_ROOT, "sensi_correct_2.sdf")
sensi_correct_path_3 = os.path.join(TEST_INPUTS_ROOT, "sensi_correct_3.sdf")
sensi_correct_path_4 = os.path.join(TEST_INPUTS_ROOT, "sensi_correct_4.sdf")

sensi_case_1 = cl.Case(sdf_path=sensi_correct_path_1)
sensi_case_2 = cl.Case(sdf_path=sensi_correct_path_2)

cov_scale_correct_path_1 = os.path.join(TEST_INPUTS_ROOT, "cov_correct_1.txt")
cov_dataf_1_scale_tuple = cl.format_scale_txt_to_dataframe(cov_scale_correct_path_1)
cov_dataf_1_scale = cov_dataf_1_scale_tuple[0]  # Extract DataFrame from tuple
cov_1_nd = cl.NDCovariances(cov_scale_correct_path_1, format="coverx_text")


class TestFunctions(unittest.TestCase):

    def assertDataframeEqual(self, a, b, msg="Both Dataframes are not equals"):
        try:
            pd.testing.assert_frame_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        # Set up any necessary variables or configurations before each test
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)
        pass

    def tearDown(self):
        # Clean up any resources or configurations after each test
        pass

    def test_formatage(self):

        # --- Test format_scale_txt_to_dataframe avec fichier vide
        try:
            cl.format_scale_txt_to_dataframe(os.path.join(TEST_INPUTS_ROOT, "fichier_vide.txt"))
        except IndexError as e:
            print("Test successfull empty file :", str(e))

        # --- Test fomattage cov data - SCALE
        if len(cov_dataf_1_scale_tuple[0]["STD"]) != 0:
            print("Test successfull for reading covariances file (type COVERX txt)")

        # --- Test fomattage cov data - COMAC
        cov_comac_correct_path_2 = os.path.join(TEST_INPUTS_ROOT, "cov_correct_2_comac")
        cov_dataf_2_comac = cl.format_comac_to_dataframe(cov_comac_correct_path_2)
        print("Test successfull for reading covariances file (type COMAC)")

        # --- Test NDCovariances object creation
        cov_obj_scale = cl.NDCovariances(input_path=cov_scale_correct_path_1, format="coverx_text")
        self.assertIsNotNone(cov_obj_scale.cov_dataf)
        self.assertEqual(cov_obj_scale.format, "coverx_text")
        self.assertIsNotNone(cov_obj_scale.iso_reac_list)
        self.assertIsInstance(cov_obj_scale.iso_reac_list, list)
        print("Test successfull for creating NDCovariances object (SCALE text)")

        cov_obj_comac = cl.NDCovariances(input_path=cov_comac_correct_path_2, format="comac")
        self.assertIsNotNone(cov_obj_comac.cov_dataf)
        self.assertEqual(cov_obj_comac.format, "comac")
        self.assertIsNotNone(cov_obj_comac.iso_reac_list)
        self.assertIsInstance(cov_obj_comac.iso_reac_list, list)
        print("Test successfull for creating NDCovariances object (COMAC)")

        try:
            cov_comac_corrupted_path_1 = os.path.join(TEST_INPUTS_ROOT, "cov_corrompu_2_comac")
            cl.format_comac_to_dataframe(cov_comac_corrupted_path_1)
        except cl.EmptyParsingError as e:
            print("Test successfull for COMAC matrix corrupted :", str(e))

        # --- Test formattage sdf avec addition des sensis pour numero d'iso et de reac egaux (typiquement H1 et H1-H2O)
        sensi_df = cl.format_sensi_to_dataframe(input_sdf_path=sensi_correct_path_3, std=True)

        self.assertAlmostEqual(sensi_df["SENSI"][3][0], 6.000000e-04, places=4)
        self.assertAlmostEqual(sensi_df["SENSI"][3][1], 6.000000e-04, places=4)
        self.assertAlmostEqual(sensi_df["SENSI_INTEGRAL"][3], 1.200000e-03, places=4)
        self.assertAlmostEqual(sensi_df["SENSI_STD"][3][0], 7.071070e-06, places=8)
        self.assertAlmostEqual(sensi_df["SENSI_STD"][3][1], 7.071070e-06, places=8)
        self.assertAlmostEqual(sensi_df["STD_INTEGRAL"][3], 1.0e-05, places=5)
        self.assertAlmostEqual(sensi_df["SENSI_INTEGRAL_ABS"][3], 1.200000e-03, places=4)

        # --- export sensi
        sensi_case_1.export_to_html(output_html_path="./sensi_test.html")

        # --- Test condense sdf
        input_path = os.path.join(TEST_INPUTS_ROOT, "sensi_correct_1.sdf")
        output_path = os.path.join(TEST_INPUTS_ROOT, "sensi_correct_1_condensed.sdf")

        cl.condense_sdf(input_sdf_path=input_path, output_ebins=[2.000000e07, 1.000000e-05], output_sdf_path=output_path)
        sensi_case_1_condensed = cl.format_sensi_to_dataframe(input_sdf_path=output_path)
        wanted_df = pd.DataFrame(
            {
                "ISO": [10000, 10000, 20000],
                "REAC": [2, 4, 2],
                "SENSI": [[0.0002], [0.0006], [0.0006]],
                "SENSI_INTEGRAL": [2.000000e-04, 6.000000e-04, 6.000000e-04],
                "STD_INTEGRAL": [7.071070e-06, 7.071070e-06, 7.071070e-06],
                "SENSI_INTEGRAL_ABS": [2.000000e-04, 6.000000e-04, 6.000000e-04],
            }
        )

        self.assertDataframeEqual(wanted_df, sensi_case_1_condensed)

        print("Test successfull for SDF file condensation")

        if os.path.exists(output_path):
            os.remove(output_path)

    def test_E(self):

        E = cl.calcul_E(case_1=sensi_correct_path_1, case_2=sensi_correct_path_2)
        self.assertAlmostEqual(E, 0.7182, places=4)

        E = cl.calcul_E(case_1=sensi_case_1, case_2=sensi_case_2)
        self.assertAlmostEqual(E, 0.7182, places=4)

        E = cl.calcul_E(case_1=sensi_correct_path_1, case_2=sensi_case_2)
        self.assertAlmostEqual(E, 0.7182, places=4)

        print("Test successfull for calculation type E")

    def test_SSR(self):

        SSR = cl.calcul_SSR(appl_case=sensi_correct_path_1, bench_case=sensi_correct_path_2)
        self.assertAlmostEqual(SSR, 0.5714, places=4)

        SSR = cl.calcul_SSR(appl_case=sensi_case_1, bench_case=sensi_case_2)
        self.assertAlmostEqual(SSR, 0.5714, places=4)

        SSR = cl.calcul_SSR(appl_case=sensi_correct_path_1, bench_case=sensi_case_2)
        self.assertAlmostEqual(SSR, 0.5714, places=4)
        print("Test successfull for calculation type SS")

    def test_G_CEA(self):

        G = cl.calcul_G_CEA(appl_case=sensi_correct_path_1, bench_case=sensi_correct_path_2)
        self.assertAlmostEqual(G, 0.5714, places=4)

        G = cl.calcul_G_CEA(appl_case=sensi_case_1, bench_case=sensi_case_2)
        self.assertAlmostEqual(G, 0.5714, places=4)

        G = cl.calcul_G_CEA(appl_case=sensi_correct_path_1, bench_case=sensi_case_2)
        self.assertAlmostEqual(G, 0.5714, places=4)

        print("Test successfull for calculation type G CEA")

    def test_G(self):

        G = cl.calcul_G(case_1=sensi_correct_path_1, case_2=sensi_correct_path_2)
        self.assertAlmostEqual(G, 0.71795, places=4)  # TODO: valeur à ajuster

        G = cl.calcul_G(case_1=sensi_case_1, case_2=sensi_case_2)
        self.assertAlmostEqual(G, 0.71795, places=4)  # TODO: valeur à ajuster

        G = cl.calcul_G(case_1=sensi_correct_path_1, case_2=sensi_case_2)
        self.assertAlmostEqual(G, 0.71795, places=4)  # TODO: valeur à ajuster

        print("Test successfull for calculation type G")

    def test_Ck(self):

        Ck = cl.calcul_Ck(case_1=sensi_correct_path_1, case_2=sensi_correct_path_2, cov_data=cov_1_nd)
        self.assertAlmostEqual(Ck, 0.7807, places=4)

        Ck = cl.calcul_Ck(case_1=sensi_case_1, case_2=sensi_case_2, cov_data=cov_1_nd)
        self.assertAlmostEqual(Ck, 0.7807, places=4)

        Ck = cl.calcul_Ck(case_1=sensi_case_1, case_2=sensi_correct_path_2, cov_data=cov_1_nd)
        self.assertAlmostEqual(Ck, 0.7807, places=4)

        # Test with NDCovariances object
        cov_obj = cl.NDCovariances(input_path=cov_scale_correct_path_1, format="coverx_text")
        Ck = cl.calcul_Ck(case_1=sensi_case_1, case_2=sensi_case_2, cov_data=cov_obj)
        self.assertAlmostEqual(Ck, 0.7807, places=4)

        print("Test successfull for calculation type Ck")

    def test_prior_unc(self):

        prior_unc = cl.calcul_uncertainty(appl_case=sensi_case_1, cov_data=cov_1_nd)
        self.assertAlmostEqual(prior_unc.value, 1.12811, places=4)

        prior_unc = cl.calcul_uncertainty(appl_case=sensi_case_1, cov_data=cov_1_nd, exclude_iso=[20000], reac_list=[2, 4])
        self.assertAlmostEqual(prior_unc.value, 8.5277e-1, places=4)

        prior_unc = cl.calcul_uncertainty(appl_case=sensi_correct_path_1, cov_data=cov_1_nd)
        self.assertAlmostEqual(prior_unc.value, 1.12811, places=4)

        # Test with NDCovariances object
        cov_obj = cl.NDCovariances(input_path=cov_scale_correct_path_1, format="coverx_text")
        prior_unc = cl.calcul_uncertainty(appl_case=sensi_case_1, cov_data=cov_obj)
        self.assertAlmostEqual(prior_unc.value, 1.12811, places=4)

        prior_unc = cl.calcul_uncertainty(appl_case=sensi_correct_path_3, cov_data=cov_1_nd)
        self.assertAlmostEqual(prior_unc.value, 1.12811, places=4)

        cl.calcul_uncertainty(appl_case=sensi_correct_path_3, cov_data=cov_1_nd)
        with open(os.path.join(cl.logs.LOG_DIR, f'CALINS_{time.strftime("%Y-%m-%d")}.log'), "r") as f:
            warn_line = f.readlines()[-4]

        self.assertEqual(warn_line.split()[2], "WARNING")

        prior_unc.export_to_html(output_html_path="./prior_unc_test.html")

        print("Test successfull for calculation type a priori uncertainty")

    def test_assimilation(self):
        # --- Test mauvaise liste isotope en entree
        try:
            assim = cl.Assimilation(
                benchmarks_list=[sensi_correct_path_1],
                appl_case=sensi_correct_path_2,
                cov_data=cov_1_nd,
                iso_reac_list=[("3", "1"), ("3", "2")],
            )
        except cl.UserInputError as e:
            print("Test successfull for wrong input iso_reac_list :", str(e))

        # --- Test calcul
        assim = cl.Assimilation(benchmarks_list=[sensi_correct_path_2], appl_case=sensi_correct_path_1, cov_data=cov_1_nd)

        self.assertAlmostEqual(assim.prior_uncertainty.value, 1.12811, places=4)
        self.assertAlmostEqual(assim.bias.value, -2.423986551e-2, places=4)
        self.assertAlmostEqual(assim.post_uncertainty.value, 1.1281, places=4)

        # --- Test calcul avec tri Ck
        assim = cl.Assimilation(
            benchmarks_list=[sensi_correct_path_2, sensi_correct_path_4],
            appl_case=sensi_correct_path_1,
            cov_data=cov_1_nd,
            Ck_threshold=0.7,
        )

        self.assertAlmostEqual(assim.prior_uncertainty.value, 1.12811, places=4)
        self.assertAlmostEqual(assim.bias.value, -2.423986551e-2, places=4)
        self.assertAlmostEqual(assim.post_uncertainty.value, 1.1281, places=4)

        assim.export_to_html(output_html_path="./assim_test.html")

        # --- Test calcul avec objet NDCovariances (option par défaut recommandée)
        assim_nd = cl.Assimilation(benchmarks_list=[sensi_correct_path_2], appl_case=sensi_correct_path_1, cov_data=cov_1_nd)

        self.assertAlmostEqual(assim_nd.prior_uncertainty.value, 1.12811, places=4)
        self.assertAlmostEqual(assim_nd.bias.value, -2.423986551e-2, places=4)
        self.assertAlmostEqual(assim_nd.post_uncertainty.value, 1.1281, places=4)

        print("Test successfull for calculation type GLLSM")

    def test_USL_gllsm(self):
        """Test GLLSM-based USL calculation."""
        # --- With 1 benchmark: N < 2 -> USL not calculable
        assim = cl.Assimilation(benchmarks_list=[sensi_correct_path_2], appl_case=sensi_correct_path_1, cov_data=cov_1_nd)

        self.assertIsInstance(assim.USL_gllsm, dict)
        expected_keys = {"USL", "calculational_margin", "bias", "std", "K", "N", "p", "q"}
        self.assertEqual(set(assim.USL_gllsm.keys()), expected_keys)
        self.assertEqual(assim.USL_gllsm["N"], 1)
        self.assertIsNone(assim.USL_gllsm["USL"])
        self.assertIsNone(assim.USL_gllsm["calculational_margin"])
        self.assertIsNone(assim.USL_gllsm["K"])

        # --- With multiple benchmarks: N >= 2 -> USL calculable
        assim2 = cl.Assimilation(
            benchmarks_list=[sensi_correct_path_2, sensi_correct_path_3, sensi_correct_path_4],
            appl_case=sensi_correct_path_1,
            cov_data=cov_1_nd,
        )

        self.assertIsNotNone(assim2.USL_gllsm["USL"])
        self.assertIsNotNone(assim2.USL_gllsm["K"])
        self.assertEqual(assim2.USL_gllsm["p"], 0.95)
        self.assertEqual(assim2.USL_gllsm["q"], 0.95)

        # Check formula consistency: USL = 1 - CM - MOS
        self.assertAlmostEqual(
            assim2.USL_gllsm["USL"],
            1.0 - assim2.USL_gllsm["calculational_margin"] - assim2.MOS,
            places=10,
        )

        # Check USL value (placeholder - à remplacer par valeur de référence)
        self.assertAlmostEqual(assim2.USL_gllsm["USL"], 0.95, places=4)
        self.assertAlmostEqual(assim2.USL_gllsm["K"], 7.6559, places=4)

        print("Test successfull for USL GLLSM method")

    def test_USL_parametric(self):
        """Test parametric USL calculation."""
        # --- With 1 benchmark: N < 2 -> USL not calculable
        assim = cl.Assimilation(benchmarks_list=[sensi_correct_path_2], appl_case=sensi_correct_path_1, cov_data=cov_1_nd)

        self.assertIsInstance(assim.USL_parametric, dict)
        expected_keys = {
            "USL", "calculational_margin", "beta", "sigma_beta", "kappa",
            "k_bar", "Delta_m", "s_k", "sigma_bar_k", "N", "p", "q",
            "shapiro_stat", "shapiro_pvalue", "normality_passed",
        }
        self.assertEqual(set(assim.USL_parametric.keys()), expected_keys)
        self.assertEqual(assim.USL_parametric["N"], 1)
        self.assertIsNone(assim.USL_parametric["USL"])
        self.assertIsNone(assim.USL_parametric["calculational_margin"])
        self.assertIsNone(assim.USL_parametric["kappa"])
        self.assertIsNone(assim.USL_parametric["sigma_beta"])
        self.assertIsNotNone(assim.USL_parametric["beta"])  # beta is still defined
        self.assertEqual(assim.USL_parametric["Delta_m"], max(0.0, assim.USL_parametric["beta"]))

        # --- With multiple benchmarks: N >= 2 -> USL calculable
        assim2 = cl.Assimilation(
            benchmarks_list=[sensi_correct_path_2, sensi_correct_path_3, sensi_correct_path_4],
            appl_case=sensi_correct_path_1,
            cov_data=cov_1_nd,
        )

        self.assertIsNotNone(assim2.USL_parametric["USL"])
        self.assertIsNotNone(assim2.USL_parametric["kappa"])
        self.assertEqual(assim2.USL_parametric["p"], 0.99)
        self.assertEqual(assim2.USL_parametric["q"], 0.99)

        # Check formula consistency: USL = 1 - CM - MOS
        self.assertAlmostEqual(
            assim2.USL_parametric["USL"],
            1.0 - assim2.USL_parametric["calculational_margin"] - assim2.MOS,
            places=10,
        )

        # Check CM formula: CM = -beta + kappa * sigma_beta + Delta_m
        self.assertAlmostEqual(
            assim2.USL_parametric["calculational_margin"],
            -assim2.USL_parametric["beta"] + assim2.USL_parametric["kappa"] * assim2.USL_parametric["sigma_beta"] + assim2.USL_parametric["Delta_m"],
            places=10,
        )

        # Check Delta_m = max(0, beta)
        self.assertEqual(assim2.USL_parametric["Delta_m"], max(0.0, assim2.USL_parametric["beta"]))

        # Check USL values (placeholders - à remplacer par valeurs de référence)
        self.assertAlmostEqual(assim2.USL_parametric["USL"], 0.90167, places=4)
        self.assertAlmostEqual(assim2.USL_parametric["beta"], 0.010050, places=6)
        self.assertAlmostEqual(assim2.USL_parametric["kappa"], 23.8956, places=4)

        # Shapiro-Wilk: with >= 3 benchmarks, normality test should produce a result
        self.assertIsNotNone(assim2.USL_parametric["shapiro_pvalue"])

        print("Test successfull for USL parametric method")

    def test_USL_parametric_normality_test(self):
        """Test that Shapiro-Wilk normality test runs with enough benchmarks."""
        assim = cl.Assimilation(
            benchmarks_list=[sensi_correct_path_2, sensi_correct_path_4, sensi_correct_path_3],
            appl_case=sensi_correct_path_1,
            cov_data=cov_1_nd,
        )

        # With >= 3 benchmarks, normality test should produce a result
        self.assertIsNotNone(assim.USL_parametric["shapiro_pvalue"])
        self.assertIsNotNone(assim.USL_parametric["shapiro_stat"])
        self.assertIsInstance(assim.USL_parametric["normality_passed"], bool)

        # Check pvalue is between 0 and 1
        self.assertGreaterEqual(assim.USL_parametric["shapiro_pvalue"], 0.0)
        self.assertLessEqual(assim.USL_parametric["shapiro_pvalue"], 1.0)

        print("Test successfull for USL parametric normality test")

    def test_USL_nonparametric(self):
        """Test nonparametric rank-order USL calculation."""
        assim = cl.Assimilation(benchmarks_list=[sensi_correct_path_2], appl_case=sensi_correct_path_1, cov_data=cov_1_nd)

        # Check that USL_nonparametric dict is populated with expected keys
        self.assertIsInstance(assim.USL_nonparametric, dict)
        expected_keys = {
            "USL", "calculational_margin", "beta", "sigma_beta",
            "min_k_tilde", "Delta_m", "sigma_worst", "CNP", "m_NP",
            "N", "p_pop", "n_sigma",
        }
        self.assertEqual(set(assim.USL_nonparametric.keys()), expected_keys)

        # Check default parameters
        self.assertEqual(assim.USL_nonparametric["p_pop"], 0.95)
        self.assertEqual(assim.USL_nonparametric["n_sigma"], 2.6)
        self.assertEqual(assim.USL_nonparametric["N"], 1)

        # With 1 benchmark: CNP = 1 - 0.95^1 = 0.05, which is <= 0.4 -> USL should be None
        self.assertAlmostEqual(assim.USL_nonparametric["CNP"], 0.05, places=4)
        self.assertIsNone(assim.USL_nonparametric["USL"])
        self.assertIsNone(assim.USL_nonparametric["calculational_margin"])
        self.assertIsNone(assim.USL_nonparametric["m_NP"])

        # Check beta = min(k_tilde) - 1
        self.assertAlmostEqual(
            assim.USL_nonparametric["beta"],
            assim.USL_nonparametric["min_k_tilde"] - 1.0,
            places=10,
        )

        # Check sigma_beta = n_sigma * sigma_worst
        self.assertAlmostEqual(
            assim.USL_nonparametric["sigma_beta"],
            assim.USL_nonparametric["n_sigma"] * assim.USL_nonparametric["sigma_worst"],
            places=10,
        )

        print("Test successfull for USL nonparametric method")

    def test_USL_nonparametric_sufficient_benchmarks(self):
        """Test nonparametric USL when enough benchmarks are provided (CNP > 0.4)."""
        # Use multiple benchmarks to get CNP > 0.4
        bench_list = [sensi_correct_path_2, sensi_correct_path_3, sensi_correct_path_4] * 5  # 15 benchmarks
        assim = cl.Assimilation(
            benchmarks_list=bench_list,
            appl_case=sensi_correct_path_1,
            cov_data=cov_1_nd,
        )

        # CNP = 1 - 0.95^N should be > 0.4 with enough benchmarks
        N_active = assim.USL_nonparametric["N"]
        expected_CNP = 1.0 - 0.95**N_active
        self.assertAlmostEqual(assim.USL_nonparametric["CNP"], expected_CNP, places=10)

        if expected_CNP > 0.4:
            # USL should be calculable
            self.assertIsNotNone(assim.USL_nonparametric["USL"])
            self.assertIsNotNone(assim.USL_nonparametric["calculational_margin"])
            self.assertIsNotNone(assim.USL_nonparametric["m_NP"])

            # Check formula consistency: USL = 1 - CM - MOS
            self.assertAlmostEqual(
                assim.USL_nonparametric["USL"],
                1.0 - assim.USL_nonparametric["calculational_margin"] - assim.MOS,
                places=10,
            )

            # Check CM = -beta + sigma_beta + Delta_m + m_NP
            self.assertAlmostEqual(
                assim.USL_nonparametric["calculational_margin"],
                -assim.USL_nonparametric["beta"] + assim.USL_nonparametric["sigma_beta"] + assim.USL_nonparametric["Delta_m"] + assim.USL_nonparametric["m_NP"],
                places=10,
            )

        print("Test successfull for USL nonparametric with sufficient benchmarks")

    def test_USL_custom_MOS(self):
        """Test that custom MOS value is correctly applied to all USL methods."""
        custom_MOS = 0.02
        assim = cl.Assimilation(
            benchmarks_list=[sensi_correct_path_2, sensi_correct_path_3, sensi_correct_path_4],
            appl_case=sensi_correct_path_1,
            cov_data=cov_1_nd,
            MOS=custom_MOS,
        )

        self.assertEqual(assim.MOS, custom_MOS)

        # GLLSM: USL = 1 - CM - MOS
        self.assertIsNotNone(assim.USL_gllsm["USL"])
        self.assertAlmostEqual(
            assim.USL_gllsm["USL"],
            1.0 - assim.USL_gllsm["calculational_margin"] - custom_MOS,
            places=10,
        )

        # Parametric: USL = 1 - CM - MOS
        self.assertIsNotNone(assim.USL_parametric["USL"])
        self.assertAlmostEqual(
            assim.USL_parametric["USL"],
            1.0 - assim.USL_parametric["calculational_margin"] - custom_MOS,
            places=10,
        )

        print("Test successfull for USL with custom MOS")

    def test_USL_parametric_custom_params(self):
        """Test parametric USL with custom p and q values."""
        assim = cl.Assimilation(
            benchmarks_list=[sensi_correct_path_2, sensi_correct_path_3, sensi_correct_path_4],
            appl_case=sensi_correct_path_1,
            cov_data=cov_1_nd,
        )

        # Store the default kappa (p=0.99, q=0.99)
        kappa_default = assim.USL_parametric["kappa"]

        # Re-run with less conservative parameters
        result = assim.calcul_USL_parametric(p=0.95, q=0.95)

        self.assertEqual(result["p"], 0.95)
        self.assertEqual(result["q"], 0.95)

        # kappa should be smaller with lower p and q (less conservative)
        self.assertLess(result["kappa"], kappa_default)

        print("Test successfull for USL parametric custom parameters")


if __name__ == "__main__":
    unittest.main()
