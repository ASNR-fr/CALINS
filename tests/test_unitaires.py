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
cov_1_nd = cl.NDCovariances(cov_scale_correct_path_1, format='coverx_text')


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
        cov_obj_scale = cl.NDCovariances(input_path=cov_scale_correct_path_1, format='coverx_text')
        self.assertIsNotNone(cov_obj_scale.cov_dataf)
        self.assertEqual(cov_obj_scale.format, 'coverx_text')
        self.assertIsNotNone(cov_obj_scale.iso_reac_list)
        self.assertIsInstance(cov_obj_scale.iso_reac_list, list)
        print("Test successfull for creating NDCovariances object (SCALE text)")
        
        cov_obj_comac = cl.NDCovariances(input_path=cov_comac_correct_path_2, format='comac')
        self.assertIsNotNone(cov_obj_comac.cov_dataf)
        self.assertEqual(cov_obj_comac.format, 'comac')
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

        SSR = cl.calcul_SSR(study_case=sensi_correct_path_1, bench_case=sensi_correct_path_2)
        self.assertAlmostEqual(SSR, 0.5714, places=4)

        SSR = cl.calcul_SSR(study_case=sensi_case_1, bench_case=sensi_case_2)
        self.assertAlmostEqual(SSR, 0.5714, places=4)

        SSR = cl.calcul_SSR(study_case=sensi_correct_path_1, bench_case=sensi_case_2)
        self.assertAlmostEqual(SSR, 0.5714, places=4)
        print("Test successfull for calculation type SS")

    def test_G(self):

        G = cl.calcul_G(study_case=sensi_correct_path_1, bench_case=sensi_correct_path_2)
        self.assertAlmostEqual(G, 0.5714, places=4)

        G = cl.calcul_G(study_case=sensi_case_1, bench_case=sensi_case_2)
        self.assertAlmostEqual(G, 0.5714, places=4)

        G = cl.calcul_G(study_case=sensi_correct_path_1, bench_case=sensi_case_2)
        self.assertAlmostEqual(G, 0.5714, places=4)

        print("Test successfull for calculation type G")

    def test_Ck(self):

        Ck = cl.calcul_Ck(case_1=sensi_correct_path_1, case_2=sensi_correct_path_2, cov_data=cov_1_nd)
        self.assertAlmostEqual(Ck, 0.7807, places=4)

        Ck = cl.calcul_Ck(case_1=sensi_case_1, case_2=sensi_case_2, cov_data=cov_1_nd)
        self.assertAlmostEqual(Ck, 0.7807, places=4)

        Ck = cl.calcul_Ck(case_1=sensi_case_1, case_2=sensi_correct_path_2, cov_data=cov_1_nd)
        self.assertAlmostEqual(Ck, 0.7807, places=4)
        
        # Test with NDCovariances object
        cov_obj = cl.NDCovariances(input_path=cov_scale_correct_path_1, format='coverx_text')
        Ck = cl.calcul_Ck(case_1=sensi_case_1, case_2=sensi_case_2, cov_data=cov_obj)
        self.assertAlmostEqual(Ck, 0.7807, places=4)

        print("Test successfull for calculation type Ck")

    def test_prior_unc(self):

        prior_unc = cl.calcul_uncertainty(study_case=sensi_case_1, cov_data=cov_1_nd)
        self.assertAlmostEqual(prior_unc.value, 1.12811, places=4)

        prior_unc = cl.calcul_uncertainty(study_case=sensi_case_1, cov_data=cov_1_nd, exclude_iso=[20000], reac_list=[2, 4])
        self.assertAlmostEqual(prior_unc.value, 8.5277e-1, places=4)

        prior_unc = cl.calcul_uncertainty(study_case=sensi_correct_path_1, cov_data=cov_1_nd)
        self.assertAlmostEqual(prior_unc.value, 1.12811, places=4)
        
        # Test with NDCovariances object
        cov_obj = cl.NDCovariances(input_path=cov_scale_correct_path_1, format='coverx_text')
        prior_unc = cl.calcul_uncertainty(study_case=sensi_case_1, cov_data=cov_obj)
        self.assertAlmostEqual(prior_unc.value, 1.12811, places=4)

        prior_unc = cl.calcul_uncertainty(study_case=sensi_correct_path_3, cov_data=cov_1_nd)
        self.assertAlmostEqual(prior_unc.value, 1.12811, places=4)

        cl.calcul_uncertainty(study_case=sensi_correct_path_3, cov_data=cov_1_nd)
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
                study_case=sensi_correct_path_2,
                cov_data=cov_1_nd,
                iso_reac_list=[("3", "1"), ("3", "2")],
            )
        except cl.UserInputError as e:
            print("Test successfull for wrong input iso_reac_list :", str(e))

        # --- Test calcul
        assim = cl.Assimilation(benchmarks_list=[sensi_correct_path_2], study_case=sensi_correct_path_1, cov_data=cov_1_nd)

        self.assertAlmostEqual(assim.prior_uncertainty.value, 1.12811, places=4)
        self.assertAlmostEqual(assim.bias.value, -2.423986551e-2, places=4)
        self.assertAlmostEqual(assim.post_uncertainty.value, 1.1281, places=4)

        # --- Test calcul avec tri Ck
        assim = cl.Assimilation(
            benchmarks_list=[sensi_correct_path_2, sensi_correct_path_4],
            study_case=sensi_correct_path_1,
            cov_data=cov_1_nd,
            Ck_threshold=0.7,
        )

        self.assertAlmostEqual(assim.prior_uncertainty.value, 1.12811, places=4)
        self.assertAlmostEqual(assim.bias.value, -2.423986551e-2, places=4)
        self.assertAlmostEqual(assim.post_uncertainty.value, 1.1281, places=4)

        assim.export_to_html(output_html_path="./assim_test.html")

        # --- Test calcul avec objet NDCovariances (option par défaut recommandée)
        assim_nd = cl.Assimilation(benchmarks_list=[sensi_correct_path_2], study_case=sensi_correct_path_1, cov_data=cov_1_nd)

        self.assertAlmostEqual(assim_nd.prior_uncertainty.value, 1.12811, places=4)
        self.assertAlmostEqual(assim_nd.bias.value, -2.423986551e-2, places=4)
        self.assertAlmostEqual(assim_nd.post_uncertainty.value, 1.1281, places=4)

        print("Test successfull for calculation type GLLSM")


if __name__ == "__main__":
    unittest.main()
