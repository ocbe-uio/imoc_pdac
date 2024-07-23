import os
from os.path import dirname
import pandas as pd
import json


class LoadDataset:

    @staticmethod
    def load_bbcsport(return_y: bool = False, return_metadata: bool = False):
        r"""
        The BBCSport dataset comprises five kinds of sports news articles (i.e., athletics, cricket, football, rugby,
        and tennis) collected from the BBC Sport website. This specific subset includes 116 samples and is
        represented by four different views.

        Samples: 116; Views: 4; Features: [1991, 2063, 2113, 2158]; Clusters: 5.

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        References
        ----------
        [paper1] D. Greene and P. Cunningham, “Practical solutions to the problem of diagonal dominance in kernel
        document clustering,” in ICML. ACM, 2006, pp. 377–384.
        [paper2] N. Rai, S. Negi, S. Chaudhury, and O. Deshmukh, “Partial multi-view clustering using graph
        regularized nmf,” in ICPR. IEEE, 2016, pp.2192–2197.
        [paper3] Jie Wen, Zheng Zhang, Lunke Fei, Bob Zhang, Yong Xu, Zhao Zhang, Jinxing Li, A Survey on Incomplete
        Multi-view Clustering, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS: SYSTEMS, 2022.
        [url] https://github.com/GPMVCDummy/GPMVC/tree/master/partialMV/PVC/recreateResults/data

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_bbcsport()
        """

        output = LoadDataset.load_dataset(dataset_name= "bbcsport", return_y = return_y, return_metadata = return_metadata)
        return output


    @staticmethod
    def load_bdgp(return_y: bool = False, return_metadata: bool = False):
        r"""
        The BDGP (Berkeley Drosophila Genome Project) is an image dataset focused on drosophila embryos,
        containing 2,500 samples featuring five distinct objects. Each sample within this dataset is characterized
        by a 1750-dimensional visual feature and a 79-dimensional textual feature.

        Samples: 2500; Views: 2; Features: [1750, 79]; Clusters: 5.

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        References
        ----------
        [paper1] X. Cai, H. Wang, H. Huang, and C. Ding, “Joint stage recognition and anatomical annotation of
        drosophila gene expression patterns,” Bioinformatics, vol. 28, no. 12, pp. i16–i24, 2012.
        [paper2] Jie Wen, Zheng Zhang, Lunke Fei, Bob Zhang, Yong Xu, Zhao Zhang, Jinxing Li, A Survey on Incomplete
        Multi-view Clustering, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS: SYSTEMS, 2022.

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_bdgp()
        """
        output = LoadDataset.load_dataset(dataset_name="bdgp", return_y=return_y, return_metadata=return_metadata)
        return output


    @staticmethod
    def load_buaa(return_y: bool = False, return_metadata: bool = False):
        r"""
        The BUAA dataset comprises facial images captured through both visual and near-infrared cameras, essentially
        representing two perspectives of the same individual. A subset of 90 visual and near infrared images from the
        first 10 classes has been selected.

        Samples: 90; Views: 2; Features: [100, 100]; Clusters: 10.

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        References
        ----------
        [paper1] D. Huang, J. Sun, and Y. Wang, “The buaa-visnir face dataset instructions,” School Comput. Sci. Eng.,
        Beihang Univ., Beijing, China, Tech. Rep. IRIP-TR-12-FR-001, 2012.
        [paper2] H. Zhao, H. Liu, and Y. Fu, “Incomplete multi-modal visual data grouping.” in IJCAI, 2016,
        pp. 2392–2398.
        [paper3] Jie Wen, Zheng Zhang, Lunke Fei, Bob Zhang, Yong Xu, Zhao Zhang, Jinxing Li, A Survey on Incomplete
        Multi-view Clustering, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS: SYSTEMS, 2022.

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_buaa()
        """
        output = LoadDataset.load_dataset(dataset_name="buaa", return_y=return_y, return_metadata=return_metadata)
        return output


    @staticmethod
    def load_caltech101(return_y: bool = False, return_metadata: bool = False):
        r"""
        The Caltech101 dataset is a widely used collection of objects, comprising 9,144 images spread across 102
        categories, which include a background category and 101 distinct objects like airplanes, ants, bass, and
        beavers. Four types of feature sets have been selected, namely Cenhist, Hog, Gist, and LBP.

        Samples: 9144; Views: 6; Features: [48, 40, 254, 1984, 512, 928]; Clusters: 102.

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        References
        ----------
        [paper1] L. Fei-Fei, R. Fergus, and P. Perona, “Learning generative visual models from few training
        examples: An incremental bayesian approach tested on 101 object categories,” in CVPR Workshop. IEEE,
        2004, pp. 178–178.
        [paper2] Li, F. Nie, H. Huang, and J. Huang, “Large-scale multi-view spectral clustering via bipartite
        graph,” in AAAI, 2015, pp. 2750–2756.
        [paper3] Jie Wen, Zheng Zhang, Lunke Fei, Bob Zhang, Yong Xu, Zhao Zhang, Jinxing Li, A Survey on Incomplete
        Multi-view Clustering, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS: SYSTEMS, 2022.
        [url] https://drive.google.com/drive/folders/1O3YmthAZGiq1ZPSdE74R7Nwos2PmnHH

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_caltech101()
        """
        output = LoadDataset.load_dataset(dataset_name= "caltech101", return_y = return_y, return_metadata = return_metadata)
        return output


    @staticmethod
    def load_digits(return_y: bool = False, return_metadata: bool = False):
        r"""
        This UCI multiple dataset showcases handwritten digit images categorized into classes 0-9, encompassing six
        distinct views. Each class comprises 200 labeled examples.

        Samples: 2000; Views: 6; Features: [6, 64, 216, 47, 76, 240]; Clusters: 10.

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        References
        ----------
        [paper1] M. van Breukelen, et al. "Handwritten digit recognition by combined classifiers", Kybernetika, 34(4):381-386,
        1998.
        [paper2] Perry, Ronan, et al. "mvlearn: Multiview Machine Learning in Python." Journal of Machine Learning Research
        22.109 (2021): 1-7.
        [paper3] Jie Wen, Zheng Zhang, Lunke Fei, Bob Zhang, Yong Xu, Zhao Zhang, Jinxing Li, A Survey on Incomplete
        Multi-view Clustering, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS: SYSTEMS, 2022.
        [url] Dheeru Dua and Casey Graff. UCI machine learning repository, 2017. URL http://archive.ics.uci.edu/ml.

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_digits()
        """
        output = LoadDataset.load_dataset(dataset_name= "digits", return_y = return_y, return_metadata = return_metadata)
        return output


    @staticmethod
    def load_metabric(return_y: bool = False, return_metadata: bool = False):
        r"""
        METABRIC offers 1904 breast cancer cases, including transcriptomic and genomic information. The labels
        indicate the cancer subtype.

        Samples: 1904; Views: 2; Features: [2000, 2000]; Clusters: 7.

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        References
        ----------
        [paper] Curtis, C., Shah, S., Chin, SF. et al. The genomic and transcriptomic architecture of 2,000 breast
         tumours reveals novel subgroups. Nature 486, 346–352 (2012). https://doi.org/10.1038/nature10983.

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_metabric()
        """
        output = LoadDataset.load_dataset(dataset_name= "metabric", return_y = return_y, return_metadata = return_metadata)
        return output


    @staticmethod
    def load_nuswide(return_y: bool = False, return_metadata: bool = False):
        r"""
        The NUSWIDE dataset, created by researchers from the National University of Singapore, is a collection of real-world
        web images. This is subset that includes 30,000 images across 31 distinct classes. Each image in this subset is
        represented by a combination of five types of low-level features, namely: color histogram; color correlogram; edge
        direction histogram; wavelet texture and block-wise color moments.

        Samples: 30000; Views: 5; Features: [64, 225, 144, 73, 128]; Clusters: 31.

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        References
        ----------
        [paper1] T.-S. Chua, J. Tang, R. Hong, H. Li, Z. Luo, and Y. Zheng, “Nus-wide: a real-world web image database from
        national university of singapore,” in ACM ICIVR, 2009, pp. 1–9.
        [paper2] Jie Wen, Zheng Zhang, Lunke Fei, Bob Zhang, Yong Xu, Zhao Zhang, Jinxing Li, A Survey on Incomplete
        Multi-view Clustering, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS: SYSTEMS, 2022.
        [url] https://drive.google.com/drive/folders/1O3YmthAZGiq1ZPSdE74R7Nwos2PmnHH

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_nuswide()
        """
        output = LoadDataset.load_dataset(dataset_name= "nuswide", return_y = return_y, return_metadata = return_metadata)
        return output


    def load_nutrimouse(return_y: bool = False, return_metadata: bool = False):
        r"""
        The nutrimouse dataset originates from a mouse nutrition study conducted by Pascal Martin at the Toxicology
        and Pharmacology Laboratory within the French National Institute for Agronomic Research. This dataset encompasses
        information from 40 mice, offering two distinct data modalities: expression levels of potentially relevant genes
        (comprising 120 numerical values) and concentrations of specific fatty acids (involving 21 numerical variables).
        Each mouse within this dataset is characterized by two labels: its genetic type, divided into two classes, and
        its diet, categorized into five classes.

        Samples: 40; Views: 2; Features: [120, 21]; Clusters: [2, 5].

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        References
        ----------
        [paper1] P. Martin, H. Guillou, F. Lasserre, S. Déjean, A. Lan, J-M. Pascussi, M. San Cristobal, P. Legrand,
        P. Besse, T. Pineau. "Novel aspects of PPARalpha-mediated regulation of lipid and xenobiotic metabolism revealed
        through a nutrigenomic study." Hepatology, 2007.
        [paper2] González I., Déjean S., Martin P.G.P and Baccini, A. (2008) CCA: "An R Package to Extend Canonical
        Correlation Analysis." Journal of Statistical Software, 23(12).
        [paper3] Perry, Ronan, et al. "mvlearn: Multiview Machine Learning in Python." Journal of Machine Learning
        Research 22.109 (2021): 1-7.

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_nutrimouse()
        """
        output = LoadDataset.load_dataset(dataset_name= "nutrimouse", return_y = return_y, return_metadata = return_metadata)
        return output


    @staticmethod
    def load_sensIT300(return_y: bool = False, return_metadata: bool = False):
        r"""
        The dataset was collected from distributed sensors within an intelligent transportation system. The dataset
        comprises 300 instances, categorized into three groups representing different real-life transportation modes.
        Each instance includes two types of data: sound recordings from sensors and vibration data, each containing
        50-dimensional characteristic attributes.

        Samples: 300; Views: 2; Features: [50, 50]; Clusters: 3.

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        References
        ----------
        [paper] Zhenjiao Liu, Zhikui Chen, Yue Li, Liang Zhao, Tao Yang, Reza Farahbakhsh, Noel Crespi, and Xiaodi
        Huang. 2023. IMC-NLT: Incomplete multi-view clustering by NMF and low-rank tensor. Expert Syst. Appl. 221, C
        (Jul 2023). https://doi.org/10.1016/j.eswa.2023.119742.
        [url] https://github.com/Liuzhenjiao123/multiview-data-sets/blob/master/sensIT300.mat

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_sensIT300()
        """
        output = LoadDataset.load_dataset(dataset_name= "sensIT300", return_y = return_y, return_metadata = return_metadata)
        return output


    @staticmethod
    def load_simulated_gm(return_y: bool = False, return_metadata: bool = False):
        r"""
        This dataset has been created using a Gaussian mixture model with the mv-learn library.

        Samples: 100; Views: 2; Features: [2, 2]; Clusters: 2.

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        References
        ----------
        [paper] Perry, Ronan, et al. "mvlearn: Multiview Machine Learning in Python." Journal of Machine Learning
        Research 22.109 (2021): 1-7.

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_simulated_gm()
        """
        output = LoadDataset.load_dataset(dataset_name= "simulated_gm", return_y = return_y, return_metadata = return_metadata)
        return output


    @staticmethod
    def load_simulated_InterSIM(return_y: bool = False, return_metadata: bool = False):
        r"""
        Three inter-related genomic datasets, methylation, gene expression and protein expression from the TCGA ovarian
        cancer study, have been created using the InterSIM tool.

        Samples: 500; Views: 3; Features: [367, 131, 160]; Clusters: 3.

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        References
        ----------
        [paper] Chalise P, Raghavan R, Fridley BL. InterSIM: Simulation tool for multiple integrative 'omic datasets'.
        Comput Methods Programs Biomed. 2016 May;128:69-74. doi: 10.1016/j.cmpb.2016.02.011. Epub 2016 Feb 27.
        PMID: 27040832; PMCID: PMC4833453.

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_simulated_InterSIM()
        """
        output = LoadDataset.load_dataset(dataset_name= "simulated_InterSIM", return_y = return_y, return_metadata = return_metadata)
        return output


    @staticmethod
    def load_simulated_netMUG(return_y: bool = False, return_metadata: bool = False):
        r"""
        The dataset consists of 1,000 synthesized samples, with two perspectives, each comprising 1,000 variables. It's
        designed to simulate intricate, interconnected datasets, such as genetic and facial data, enabling the
        exploration of complex relationships.

        Samples: 1000; Views: 2; Features: [2000, 2000]; Clusters: 3.

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        References
        ----------
        [paper] Li Z, Melograna F, Hoskens H, Duroux D, Marazita ML, Walsh S, Weinberg SM, Shriver MD, Müller-Myhsok B,
        Claes P, Van Steen K. netMUG: a novel network-guided multi-view clustering workflow for dissecting genetic and
        facial heterogeneity. Front Genet. 2023 Dec 6;14:1286800. doi: 10.3389/fgene.2023.1286800. PMID: 38125750;
        PMCID: PMC10731261.

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_simulated_netMUG()
        """
        output = LoadDataset.load_dataset(dataset_name= "simulated_netMUG", return_y = return_y, return_metadata = return_metadata)
        return output


    @staticmethod
    def load_statlog(return_y: bool = False, return_metadata: bool = False):
        r"""
        The image segmentation dataset was randomly chosen from a database encompassing images from seven distinct
        categories. These images underwent manual segmentation to assign classifications to individual pixels. Compiled
        by the vision group at the University of Massachusetts, the dataset comprises 2,310 instances, each associated
        with categories represented across two perspectives. One view has a characteristic dimension of 9, while the
        other view has a characteristic dimension of 10.

        Samples: 2310; Views: 2; Features: [9, 10]; Clusters: 7.

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        References
        ----------
        [paper] Zhenjiao Liu, Zhikui Chen, Yue Li, Liang Zhao, Tao Yang, Reza Farahbakhsh, Noel Crespi, and Xiaodi
        Huang. 2023. IMC-NLT: Incomplete multi-view clustering by NMF and low-rank tensor. Expert Syst. Appl. 221, C
        (Jul 2023). https://doi.org/10.1016/j.eswa.2023.119742.
        [url] https://github.com/Liuzhenjiao123/multiview-data-sets/tree/master

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_statlog()
        """
        output = LoadDataset.load_dataset(dataset_name= "statlog", return_y = return_y, return_metadata = return_metadata)
        return output


    @staticmethod
    def load_tcga(return_y: bool = False, return_metadata: bool = False):
        r"""
        This dataset is composed of ten cancer types multi-omics data from The Cancer Genome Atlas (TCGA). This is a subset
        composed of four kinds of data: mRNA, miRNA, DNA-methylation and proteomics. Two possible targets are provided:
        origin tissue (10 labels) and survival data (numerical values).

        Samples: 2437; Views: 4; Features: [215, 2000, 131, 1739]; Clusters: 10.

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        References
        ----------
        [paper] Hoadley, Katherine & Yau, Christina & Wolf, Denise & Cherniack, Andrew & Tamborero, David & Ng, Sam &
        Leiserson, Mark & Niu, Shubin & Mclellan, Michael & Uzunangelov, Vladislav & Zhang, Jiashan & Kandoth, Cyriac &
        Akbani, Rehan & Shen, Hui & Omberg, Larsson & Chu, Andy & Margolin, Adam & van 't Veer, Laura & López-Bigas, Nuria
        & Zou, Lihua. (2014). Multiplatform Analysis of 12 Cancer Types Reveals Molecular Classification within and across
        Tissues of Origin. Cell. 158. 10.1016/j.cell.2014.06.049.
        [url] https://www.synapse.org

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_tcga()
        """
        output = LoadDataset.load_dataset(dataset_name= "tcga", return_y = return_y, return_metadata = return_metadata)
        return output


    @staticmethod
    def load_webkb(return_y: bool = False, return_metadata: bool = False):
        r"""
        The dataset consists of course and non-course documents, each offering two representations: the textual content
        of the webpage and the anchor text linked to other webpages referencing it. With regard to the page
        representation, 3,000 features were chosen, while 1,840 features were generated for linked representations.

        Samples: 1051; Views: 2; Features: [3000, 1840]; Clusters: 2.

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        References
        ----------
        [paper] Zhenjiao Liu, Zhikui Chen, Yue Li, Liang Zhao, Tao Yang, Reza Farahbakhsh, Noel Crespi, and Xiaodi
        Huang. 2023. IMC-NLT: Incomplete multi-view clustering by NMF and low-rank tensor. Expert Syst. Appl. 221, C
        (Jul 2023). https://doi.org/10.1016/j.eswa.2023.119742.
        [url] https://lig-membres.imag.fr/grimal/data.html

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_webkb()
        """
        output = LoadDataset.load_dataset(dataset_name= "webkb", return_y = return_y, return_metadata = return_metadata)
        return output


    @staticmethod
    def load_wisconsin(return_y: bool = False, return_metadata: bool = False):
        r"""
        The dataset comprises webpages extracted from the University of Wisconsin website, categorized into five types:
        student, project, course, staff, and faculty. Each webpage offers two perspectives: the content view and the
        reference view. In the content view, every webpage contains 1,703 words. The reference view delineates the
        reference relationships between a page and other pages within the dataset.

        Samples: 265; Views: 2; Features: [265, 1703]; Clusters: 5.

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        References
        ----------
        [paper] Zhenjiao Liu, Zhikui Chen, Yue Li, Liang Zhao, Tao Yang, Reza Farahbakhsh, Noel Crespi, and Xiaodi
        Huang. 2023. IMC-NLT: Incomplete multi-view clustering by NMF and low-rank tensor. Expert Syst. Appl. 221, C
        (Jul 2023). https://doi.org/10.1016/j.eswa.2023.119742.
        [url] https://lig-membres.imag.fr/grimal/data.html

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_wisconsin()
        """
        output = LoadDataset.load_dataset(dataset_name= "wisconsin", return_y = return_y, return_metadata = return_metadata)
        return output


    @staticmethod
    def load_dataset(dataset_name: str, return_y: bool = False, return_metadata: bool = False):
        r"""
        Load a multi-view dataset.

        Parameters
        ----------
        dataset_name: str
            Name of the dataset. It must be one of: "bbcsport", "bdgp", "buaa", "caltech101", "digits", "metabric",
            "nuswide", "nutrimouse", "simulated_gm", "simulated_InterSIM", "simulated_netMUG", "tcga".
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_dataset(dataset_name = 'tcga')
        """
        module_path = dirname(__file__)
        data_path = os.path.join(module_path, "data", dataset_name)
        data_files = [filename for filename in os.listdir(data_path)]
        data_files = sorted(data_files)
        data_files = [os.path.join(data_path, filename) for filename in data_files if dataset_name in filename and not filename.endswith("y.csv")]
        Xs = [pd.read_csv(filename) for filename in data_files]
        output = (Xs,)
        if return_y:
            y = pd.read_csv(os.path.join(data_path, f"{dataset_name}_y.csv"))
            y = y.loc[Xs[0].index]
            if y.shape[1] > 1:
                y = y.squeeze()
            output = output + (y,)
        if return_metadata:
            metadata_filename = os.path.join(data_path, "metadata.json")
            if os.path.isfile(metadata_filename):
                with open(metadata_filename) as json_file:
                    metadata = json.load(json_file)
                output = output + (metadata,)
            else:
                output = output + (None,)
        if len(output) == 1:
            output = output[0]
        return output
