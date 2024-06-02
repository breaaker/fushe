from PyQt6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QComboBox, QScrollArea)
from PyQt6.QtGui import QPixmap
import search_params as sep
import numpy as np

rho_air = 1.293e-3
rho_soft_tissue = 1

class Page_cal_out_beta(QWidget):
    ##############################
    # 外照射
    ##############################
    def __init__(self):
        super().__init__()
        self.initUI()
        
    def initUI(self):
        ##############################
        # beta
        ##############################

        layout = QHBoxLayout()
        self.setLayout(layout)
        
        ##已知注量率
        
        beta_energy = QLineEdit(self)
        beta_energy.setPlaceholderText("请输入电子的能量,单位:MeV,0.01-10MeV")
        beta_phi = QLineEdit(self)
        beta_phi.setPlaceholderText("请输入电子的注量率,单位:/(m^2s)")
        button = QPushButton("Calculate", self)
        beta_dose = QLabel()
        beta_dose.setText("剂量率: ")
        button.clicked.connect(lambda: self.cal_beta(beta_energy, beta_phi, beta_dose))

        beta_layout = QVBoxLayout()
        beta_layout.addWidget(QLabel("############################################"))
        beta_layout.addWidget(QLabel("已知注量率"))
        beta_layout.addWidget(beta_energy)
        beta_layout.addWidget(beta_phi)
        beta_layout.addWidget(beta_dose)
        beta_layout.addWidget(button)
        
        ##洛文格公式
        beta_source = QLineEdit(self)
        beta_source.setPlaceholderText("请输入源的种类,例如:Sr-90")
        beta_A = QLineEdit(self)
        beta_A.setPlaceholderText("请输入源的活度,单位:Bq")
        beta_material = QComboBox(self)
        beta_material.addItem("空气")
        beta_material.addItem("软组织")
        beta_distance = QLineEdit(self)
        beta_distance.setPlaceholderText("请输入距离,单位:cm")
        beta_dose_rate = QLabel()
        beta_dose_rate.setText("剂量率: ")
        button_2 = QPushButton("Calculate", self)
        button_2.clicked.connect(lambda: self.cal_beta_lovinger(beta_source, beta_A, beta_material, beta_distance, beta_dose_rate))

        beta_layout.addWidget(QLabel("############################################"))
        beta_layout.addWidget(QLabel("beta点源,洛文格公式"))
        beta_layout.addWidget(beta_source)
        beta_layout.addWidget(beta_A)
        beta_layout.addWidget(beta_material)
        beta_layout.addWidget(beta_distance)
        beta_layout.addWidget(beta_dose_rate)
        beta_layout.addWidget(button_2)

        ##积分洛文格公式
        beta_source_2 = QLineEdit(self)
        beta_source_2.setPlaceholderText("请输入源的种类,例如:Sr-90")
        beta_A_2 = QLineEdit(self)
        beta_A_2.setPlaceholderText("请输入源的活度面密度,单位:Bq/cm^2")
        beta_material_2 = QComboBox(self)
        beta_material_2.addItem("空气")
        beta_material_2.addItem("软组织")
        beta_distance_2 = QLineEdit(self)
        beta_distance_2.setPlaceholderText("请输入距离,单位:cm")
        beta_dose_rate_2 = QLabel()
        beta_dose_rate_2.setText("剂量率: ")
        button_3 = QPushButton("Calculate", self)
        button_3.clicked.connect(lambda: self.cal_beta_lovinger(beta_source_2, beta_A_2, beta_material_2, beta_distance_2, beta_dose_rate_2, True))

        beta_layout_2 = QVBoxLayout()
        beta_layout_2.addWidget(QLabel("############################################"))
        beta_layout_2.addWidget(QLabel("beta面源,积分洛文格公式"))
        beta_layout_2.addWidget(beta_source_2)
        beta_layout_2.addWidget(beta_A_2)
        beta_layout_2.addWidget(beta_material_2)
        beta_layout_2.addWidget(beta_distance_2)
        beta_layout_2.addWidget(beta_dose_rate_2)
        beta_layout_2.addWidget(button_3)

        layout.addLayout(beta_layout)
        layout.addLayout(beta_layout_2)

    def cal_beta(self, beta_energy, beta_phi, dose):
        '''
        该函数用于计算beta辐射的剂量
        '''
        phi_1 = float(beta_phi.text())
        phi_0 = QLabel()
        phi_0 = float(sep.Page_search_beta().search_phi(beta_energy, phi_0, True))
        re = 2.5e-3 * phi_1 / phi_0
        dose.setText("剂量率: " + str(re) + "mGy/h")
    
    def cal_beta_lovinger(self, beta_source, beta_A, beta_material, beta_distance, beta_dose_rate, surface = False):
        '''
        该函数用于计算beta辐射的剂量
        '''
        A = float(beta_A.text())
        over = 1
        if beta_source.text() == "Sr-90":
            over = 1.17
        
        information = sep.Page_search_beta().search_source(beta_source, QLabel(), True)
        E_max = float(information[1])
        if E_max < 0.167 or E_max > 2.24:
            beta_dose_rate.setText("超出洛文格公式适用范围")
            return
        E_mean = float(information[2])

        if beta_material.currentText() == "空气":
            rho = rho_air
            c = 3.11*np.e**(-0.56*E_max)
            v = 16/(E_max-0.036)**1.4 *(2-over)
        else:
            rho = rho_soft_tissue
            if E_max < 0.5:
                c = 2
            elif E_max < 1.5:
                c = 1.5
            else:
                c = 1
            v = 18.6/(E_max-0.036)**1.37 *(2-over)
        alpha = 1/(3*c**2 - np.e*(c**2-1))
        K = 4.59e-5 * rho**2 * v**3 * E_mean * alpha
        r = float(beta_distance.text())*rho
        if surface:
            D = 2.89e-4 * v * E_mean * alpha *A*(c*(1 + np.log(c/(v*r)) - np.e**(1-v*r/c)) + np.e**(1-v*r))
        else:
            if v*r/c >= 1:
                D = K * A * np.e**(1-v*r) / (v*r)
            else:
                D = K*A/(v*r)**2 * (c*(1-v*r/c * np.e**(1-v*r/c)) + v*r*np.e**(1-v*r))
        beta_dose_rate.setText("剂量率: " + str(D) + "mGy/h")

with open("data/gamma_data/point.txt", "r") as f:
    data = f.readlines()
gamma_point = {}
for i in range(len(data)):
    gamma_data = data[i].strip().split("\t")
    gamma_point[gamma_data[0]] = [float(gamma_data[1]), float(gamma_data[2])]

class Page_cal_out_gamma(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
    def initUI(self):
        ##############################
        # gamma
        ##############################

        ## 点源剂量率
        layout = QHBoxLayout()
        self.setLayout(layout)
        
        gamma_source = QLineEdit(self)
        gamma_source.setPlaceholderText("请输入源的种类,例如:Co-60")
        gamma_A = QLineEdit(self)
        gamma_A.setPlaceholderText("请输入源的活度,单位:Bq")
        gamma_distance = QLineEdit(self)
        gamma_distance.setPlaceholderText("请输入距离,单位:m")
        gamma_X = QLabel()
        gamma_X.setText("照射量率: ")
        gamma_dose_rate = QLabel()
        gamma_dose_rate.setText("剂量率: ")
        gamma_K = QLabel()
        gamma_K.setText("比释动能率: ")
        button = QPushButton("Calculate", self)
        button.clicked.connect(lambda: self.cal_gamma_point(gamma_source, gamma_A, gamma_distance, gamma_X, gamma_dose_rate, gamma_K))

        gamma_layout = QVBoxLayout()
        gamma_layout.addWidget(QLabel("Gamma part"))
        gamma_layout.addWidget(QLabel("############################################"))
        gamma_layout.addWidget(QLabel("点源剂量率"))
        gamma_layout.addWidget(gamma_source)
        gamma_layout.addWidget(gamma_A)
        gamma_layout.addWidget(gamma_distance)
        gamma_layout.addWidget(gamma_X)
        gamma_layout.addWidget(gamma_dose_rate)
        gamma_layout.addWidget(gamma_K)
        gamma_layout.addWidget(button)

        gamma_source_2 = QLineEdit(self)
        gamma_source_2.setPlaceholderText("请输入源的种类,例如:Co-60")
        gamma_A_2 = QLineEdit(self)
        gamma_A_2.setPlaceholderText("请输入源的线活度,单位:Bq/m")
        start = QLineEdit(self)
        start.setPlaceholderText("请输入起始坐标, 例如(1,0,0), 单位:m")
        end = QLineEdit(self)
        end.setPlaceholderText("请输入终止坐标, 例如(0,1,0), 单位:m")
        point = QLineEdit(self)
        point.setPlaceholderText("请输入测量点的坐标, 例如(0,0,0), 单位:m")
        gamma_X_2 = QLabel()
        gamma_X_2.setText("照射量率: ")
        gamma_dose_rate_2 = QLabel()
        gamma_dose_rate_2.setText("剂量率: ")
        gamma_K_2 = QLabel()
        gamma_K_2.setText("比释动能率: ")
        button_2 = QPushButton("Calculate", self)
        button_2.clicked.connect(lambda: self.cal_gamma_line(gamma_source_2, gamma_A_2, start, end, point, gamma_X_2, gamma_dose_rate_2, gamma_K_2))

        gamma_layout.addWidget(QLabel("############################################"))
        gamma_layout.addWidget(QLabel("线源剂量率"))
        gamma_layout.addWidget(gamma_source_2)
        gamma_layout.addWidget(gamma_A_2)
        x_layout = QHBoxLayout()
        x_layout.addWidget(start)
        x_layout.addWidget(end)
        gamma_layout.addLayout(x_layout)
        gamma_layout.addWidget(point)
        gamma_layout.addWidget(gamma_X_2)
        gamma_layout.addWidget(gamma_dose_rate_2)
        gamma_layout.addWidget(gamma_K_2)
        gamma_layout.addWidget(button_2)

        gamma_source_3 = QLineEdit(self)
        gamma_source_3.setPlaceholderText("请输入源的种类,例如:Co-60")
        gamma_A_3 = QLineEdit(self)
        gamma_A_3.setPlaceholderText("请输入源的面活度,单位:Bq/m^2")
        r = QLineEdit(self)
        r.setPlaceholderText("请输入面源的半径,单位:m")
        h = QLineEdit(self)
        h.setPlaceholderText("请输入测量点到面源的距离,单位:m")
        gamma_X_3 = QLabel()
        gamma_X_3.setText("照射量率: ")
        gamma_dose_rate_3 = QLabel()
        gamma_dose_rate_3.setText("剂量率: ")
        gamma_K_3 = QLabel()
        gamma_K_3.setText("比释动能率: ")
        button_3 = QPushButton("Calculate", self)
        button_3.clicked.connect(lambda: self.cal_gamma_surface(gamma_source_3, gamma_A_3, r, h, gamma_X_3, gamma_dose_rate_3, gamma_K_3))

        gamma_layout_2 = QVBoxLayout()
        gamma_layout_2.addWidget(QLabel("############################################"))
        gamma_layout_2.addWidget(QLabel("面源剂量率"))
        gamma_layout_2.addWidget(gamma_source_3)
        gamma_layout_2.addWidget(gamma_A_3)
        x_layout_2 = QHBoxLayout()
        x_layout_2.addWidget(r)
        x_layout_2.addWidget(h)
        gamma_layout_2.addLayout(x_layout_2)
        gamma_layout_2.addWidget(gamma_X_3)
        gamma_layout_2.addWidget(gamma_dose_rate_3)
        gamma_layout_2.addWidget(gamma_K_3)
        gamma_layout_2.addWidget(button_3)

        layout.addLayout(gamma_layout)
        layout.addLayout(gamma_layout_2)

    def cal_gamma_point(self, gamma_source, gamma_A, gamma_distance, gamma_X, gamma_dose_rate, gamma_K):
        '''
        该函数用于计算gamma辐射的剂量
        '''
        source = gamma_source.text()
        if source not in gamma_point:
            gamma_X.setText("未知的源")
            gamma_dose_rate.setText("未知的源")
            gamma_K.setText("未知的源")
            return
        Gamma = gamma_point[source][0]
        Gamma_K = gamma_point[source][1]
        A = float(gamma_A.text())
        r = float(gamma_distance.text())
        X = A * Gamma / r**2
        K = A * Gamma_K / r**2
        gamma_X.setText("照射量率: " + str(X) + "C/kg/s")
        gamma_dose_rate.setText("剂量率: " + str(X/2.58e-4*8.69e-3*3.6e3) + "Gy/h")        
        gamma_K.setText("比释动能率: " + str(K) + "Gy/s")

    def cal_gamma_line(self, gamma_source, gamma_A, start, end, point, gamma_X, gamma_dose_rate, gamma_K):
        '''
        该函数用于计算gamma辐射的剂量
        '''
        source = gamma_source.text()
        if source not in gamma_point:
            gamma_X.setText("未知的源")
            gamma_dose_rate.setText("未知的源")
            gamma_K.setText("未知的源")
            return
        Gamma = gamma_point[source][0]
        Gamma_K = gamma_point[source][1]
        A = float(gamma_A.text())

        #计算theta和h
        start = np.array(eval(start.text()))
        end = np.array(eval(end.text()))
        point = np.array(eval(point.text()))
        theta = np.arccos(np.dot(start-point, end-point)/(np.linalg.norm(start-point)*np.linalg.norm(end-point)))
        l1 = np.dot(start-point, start-end)/np.linalg.norm(start-end)
        h = np.sqrt(np.linalg.norm(start-point)**2 - l1**2)

        X = A * Gamma * theta / h
        K = A * Gamma_K * theta / h
        gamma_X.setText("照射量率: " + str(X) + "C/kg/s")
        gamma_dose_rate.setText("剂量率: " + str(X/2.58e-4*8.69e-3*3.6e3) + "Gy/h")
        gamma_K.setText("比释动能率: " + str(K) + "Gy/s")

    def cal_gamma_surface(self, gamma_source, gamma_A, r, h, gamma_X, gamma_dose_rate, gamma_K):
        '''
        该函数用于计算gamma辐射的剂量
        '''
        source = gamma_source.text()
        if source not in gamma_point:
            gamma_X.setText("未知的源")
            gamma_dose_rate.setText("未知的源")
            gamma_K.setText("未知的源")
            return
        Gamma = gamma_point[source][0]
        Gamma_K = gamma_point[source][1]
        A = float(gamma_A.text())
        r = float(r.text())
        h = float(h.text())
        X = Gamma * A * np.pi * np.log(r**2/h**2 + 1)
        K = Gamma_K * A * np.pi * np.log(r**2/h**2 + 1)
        gamma_X.setText("照射量率: " + str(X) + "C/kg/s")
        gamma_dose_rate.setText("剂量率: " + str(X/2.58e-4*8.69e-3*3.6e3) + "Gy/h")
        gamma_K.setText("比释动能率: " + str(K) + "Gy/s")

class Page_cal_out_neutron(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
    def initUI(self):
        ##############################
        # neutron
        ##############################

        ## 解析法算剂量
        layout = QHBoxLayout()
        self.setLayout(layout)

        material = QComboBox(self)
        material.addItem("水")
        material.addItem("干燥空气")
        material.addItem("标准人")
        material.addItem("肌肉(ICRU)")
        material.addItem("骨(股骨)")
        material.addItem("近似组织")
        material.addItem("尼龙")
        material.addItem("有机玻璃")

        energy = QLineEdit(self)
        energy.setPlaceholderText("请输入中子能量, 单位MeV, 范围为0.01e-4到29MeV")

        neutron_phi = QLineEdit(self)
        neutron_phi.setPlaceholderText("请输入中子注量, 单位/(cm^2)")
        neutron_K = QLabel()
        neutron_K.setText("K: ")
        button = QPushButton("Calculate", self)
        button.clicked.connect(lambda: self.cal_neutron(material, energy, neutron_phi, neutron_K))

        neutron_layout = QVBoxLayout()
        neutron_layout.addWidget(QLabel("Neutron part"))
        neutron_layout.addWidget(QLabel("############################################"))
        neutron_layout.addWidget(QLabel("解析法算剂量"))
        neutron_layout.addWidget(material)
        neutron_layout.addWidget(energy)
        neutron_layout.addWidget(neutron_phi)
        neutron_layout.addWidget(neutron_K)
        neutron_layout.addWidget(button)

        ## 经验法算剂量
        n_source = QComboBox(self)
        n_source.addItem("Po-Be")
        n_source.addItem("Ra-Be")
        n_source.addItem("Pu9-Be")
        n_source.addItem("Am-Be")

        n_A = QLineEdit(self)
        n_A.setPlaceholderText("请输入源的活度, 单位:Bq")
        n_distance = QLineEdit(self)
        n_distance.setPlaceholderText("请输入距离, 单位:m")
        n_H = QLabel()
        n_H.setText("中子剂量当量率: ")

        button_2 = QPushButton("Calculate", self)
        button_2.clicked.connect(lambda: self.cal_neutron_empirical(n_source, n_A, n_distance, n_H))

        neutron_layout.addWidget(QLabel("############################################"))
        neutron_layout.addWidget(QLabel("经验法算剂量"))
        neutron_layout.addWidget(n_source)
        neutron_layout.addWidget(n_A)
        neutron_layout.addWidget(n_distance)
        neutron_layout.addWidget(n_H)
        neutron_layout.addWidget(button_2)

        ## 中子屏蔽
        n_material = QComboBox(self)
        n_material.addItem("水")
        n_material.addItem("石蜡")
        n_material.addItem("聚乙烯")
        n_material.addItem("铅")
        n_material.addItem("铁")

        n_source_2 = QComboBox(self)
        n_source_2.addItem("Na-Be")
        n_source_2.addItem("Sb-Be")
        n_source_2.addItem("Po-Be")
        n_source_2.addItem("Ra-Be")
        n_source_2.addItem("Pu8-Be")
        n_source_2.addItem("Pu9-Be")
        n_source_2.addItem("Am-Be")

        n_A_2 = QLineEdit(self)
        n_A_2.setPlaceholderText("请输入源的活度, 单位:Bq")
        n_distance_2 = QLineEdit(self)
        n_distance_2.setPlaceholderText("请输入距离, 单位:cm")
        n_phi = QLineEdit(self)
        n_phi.setPlaceholderText("请输入容许的中子注量率, 单位/(cm^2)")
        q = QLineEdit(self)
        q.setPlaceholderText("请输入居留因子")
        n_thick = QLabel()
        n_thick.setText("所需材料的厚度: ")
        button_3 = QPushButton("Calculate", self)
        button_3.clicked.connect(lambda: self.cal_neutron_shielding(n_material, n_source_2, n_A_2, n_distance_2, n_phi, q, n_thick))

        neutron_layout_2 = QVBoxLayout()
        neutron_layout_2.addWidget(QLabel("############################################"))
        neutron_layout_2.addWidget(QLabel("中子屏蔽"))
        neutron_layout_2.addWidget(n_material)
        neutron_layout_2.addWidget(n_source_2)
        neutron_layout_2.addWidget(n_A_2)
        neutron_layout_2.addWidget(n_distance_2)
        neutron_layout_2.addWidget(n_phi)
        neutron_layout_2.addWidget(q)
        neutron_layout_2.addWidget(n_thick)
        neutron_layout_2.addWidget(button_3)

        layout.addLayout(neutron_layout)
        layout.addLayout(neutron_layout_2)
        
    def cal_neutron(self, material, energy, neutron_phi, neutron_K):
        '''
        该函数用于计算neutron辐射的剂量
        '''
        f_k = sep.Page_search_neutron().search_f_k(material.currentText(), energy.text(), QLabel(), True)
        K = float(neutron_phi.text()) * f_k
        neutron_K.setText("K: " + str(K) + "Gy")

    def cal_neutron_empirical(self, n_source, n_A, n_distance, n_H):
        '''
        该函数用于计算neutron辐射的剂量
        '''
        data = sep.Page_search_neutron().search_n_source(n_source, QLabel(), True)
        Y = float(data[2])
        f = float(data[3])
        A = float(n_A.text())
        r = float(n_distance.text())
        H = A * Y / (4 * np.pi * r**2) * f
        n_H.setText("中子剂量当量率: " + str(H) + "Sv/s")

    def cal_neutron_shielding(self, n_material, n_source, n_A, n_distance, n_phi, q, n_thick):
        material = n_material.currentText()
        if material == "铅":
            B=3.5
        elif material == "铁":
            B=2.6
        else:
            B=5
        data = sep.Page_search_neutron().search_n_source(n_source, QLabel(), True)
        Y = float(data[2])
        sigma = sep.Page_search_neutron().search_n_sigma(material, QLabel(), True)
        A = float(n_A.text())
        r = float(n_distance.text())
        phi = float(n_phi.text())
        q = float(q.text())
        thick = 1/sigma * np.log(A*Y*B*q/(4*np.pi*r**2*phi))
        n_thick.setText("所需材料的厚度: " + str(thick) + "cm")

from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.font_manager import FontProperties
import scipy.special._cdflib

font = FontProperties(family="SimSun", size=14)
plt.rcParams['font.sans-serif'] = ['SimSun']
plt.rcParams['axes.unicode_minus'] = False

class Page_cal_in(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        
        ##############################
        # 内照射部分
        ##############################

        ##设置隔室链, 并进行可视化
        layout = QHBoxLayout()
        self.setLayout(layout)
        
        self.Ts = []
        self.nodes = []
        self.edges = np.zeros((0, 0))
        self.As = np.zeros(0)

        ##增加节点
        G = nx.DiGraph()
        set_nuc_T = QLineEdit(self)
        set_nuc_T.setPlaceholderText("请输入核素的半衰期, 单位:s")

        button_nuc = QPushButton("设置", self)
        button_nuc.clicked.connect(lambda: self.set_T(set_nuc_T))

        add_node = QLineEdit(self)
        add_node.setPlaceholderText("请输入要添加的隔室名称, 例如:胃")
        set_T = QLineEdit(self)
        set_T.setPlaceholderText("请输入核素在该节点的生物半排期, 单位:s")

        state_node = QLabel()
        state_node.setText("状态: ")

        start = QHBoxLayout()
        start.addWidget(QLabel("添加起始的隔室："))
        add_edge_start = QComboBox(self)
        start.addWidget(add_edge_start)
        end = QHBoxLayout()
        end.addWidget(QLabel("添加终止的隔室："))
        add_edge_end = QComboBox(self)
        end.addWidget(add_edge_end)

        chain_pic = QScrollArea(self)

        button = QPushButton("添加", self)
        button.clicked.connect(lambda: self.add_nodes(G, add_node, set_T, add_edge_start, add_edge_end, chain_pic, state_node, all_node))
        
        in_layout = QVBoxLayout()
        in_layout.addWidget(QLabel("内照射部分"))
        in_layout.addWidget(QLabel("############################################"))
        x_layout_1 = QHBoxLayout()
        x_layout_1.addWidget(set_nuc_T)
        x_layout_1.addWidget(button_nuc)
        in_layout.addLayout(x_layout_1)

        x_layout_2 = QHBoxLayout()
        x_layout_2.addWidget(add_node)
        x_layout_2.addWidget(set_T)
        in_layout.addLayout(x_layout_2)
        
        in_layout.addWidget(state_node)
        in_layout.addWidget(button)

        ##增加边
        ratio = QLineEdit(self)
        ratio.setPlaceholderText("请输入该核素从起始隔室流向终止隔室的份额, 例如:0.5")
        state = QLabel()
        state.setText("状态: ")
        button_2 = QPushButton("添加", self)
        button_2.clicked.connect(lambda: self.add_edge(G, add_edge_start, add_edge_end, ratio, chain_pic, state))

        in_layout.addLayout(start)
        in_layout.addLayout(end)
        in_layout.addWidget(ratio)
        in_layout.addWidget(state)
        in_layout.addWidget(button_2)
        in_layout.addWidget(chain_pic)

        ##计算剂量
        all_node = QComboBox(self)

        origin_A = QLineEdit(self)
        origin_A.setPlaceholderText("请输入该隔室的初始活度, 单位:Bq, 不设置的话默认为0")
        state_A = QLabel()
        state_A.setText("状态: ")
        button_3 = QPushButton("添加", self)
        button_3.clicked.connect(lambda: self.set_A(all_node, origin_A, state_A))

        in_layout_2 = QVBoxLayout()
        in_layout_2.addWidget(all_node)
        in_layout_2.addWidget(origin_A)
        in_layout_2.addWidget(state_A)
        in_layout_2.addWidget(button_3)

        time = QLineEdit(self)
        time.setPlaceholderText("请输入想要计算的时长, 单位:s")
        A_pic = QScrollArea(self)
        button_4 = QPushButton("计算", self)
        button_4.clicked.connect(lambda: self.cal_in(time, A_pic))

        in_layout_2.addWidget(time)
        in_layout_2.addWidget(button_4)
        in_layout_2.addWidget(A_pic)

        layout.addLayout(in_layout)
        layout.addLayout(in_layout_2)

    def set_T(self, set_nuc_T):
        '''
        该函数用于设置核素的半衰期
        '''
        self.T = eval(set_nuc_T.text())

    def add_nodes(self, G, add_node, set_T, add_edge_start, add_edge_end, chain_pic, state, all_node):
        '''
        该函数用于添加隔室
        '''
        node = add_node.text()
        if node in self.nodes:
            state.setText("状态: 添加失败，该隔室已存在")
            return
        G.add_node(node)
        add_edge_start.addItem(node)
        add_edge_end.addItem(node)
        all_node.addItem(node)

        self.nodes.append(node)
        self.Ts.append(eval(set_T.text()))
        a = self.edges
        a = np.vstack((a, np.zeros(a.shape[0])))
        a = np.hstack((a, np.zeros((a.shape[0], 1))))
        self.edges = a
        self.As = np.append(self.As, 0)

        plt.figure(figsize=(3, 2))
        plt.clf()
        nx.draw(G, with_labels=True, node_size=1000, node_color ="skyblue", font_size=15, font_family='SimSun', arrows=True)
        plt.savefig("chain.png")
        pixmap = QPixmap("chain.png")
        pic = QLabel(self)
        pic.setPixmap(pixmap)
        chain_pic.setWidget(pic)
        state.setText("状态: 添加成功")

    def add_edge(self, G, add_edge_start, add_edge_end, ratio, chain_pic, state):
        '''
        该函数用于添加边
        '''
        start = add_edge_start.currentText()
        end = add_edge_end.currentText()
        if start == end:
            state.setText("状态: 添加失败, 起始隔室和终止隔室不能相同")
            return
        if self.edges[self.nodes.index(end)][self.nodes.index(start)] != 0:
            state.setText("状态: 添加失败, 该边已存在")
            return
        ratio = eval(ratio.text())
        G.add_edge(start, end)
        start_index = self.nodes.index(start)
        end_index = self.nodes.index(end)
        self.edges[end_index][start_index] = ratio
        plt.figure(figsize=(3, 2))
        plt.clf()
        nx.draw(G, with_labels=True, node_size=1000,node_color ="skyblue", font_size=15, font_family='SimSun', arrows=True)
        plt.savefig("chain.png")
        pixmap = QPixmap("chain.png")
        pic = QLabel(self)
        pic.setPixmap(pixmap)
        chain_pic.setWidget(pic)
        state.setText("状态: 添加成功")

    def set_A(self, node, origin_A, state_A):
        '''
        该函数用于设置初始活度
        '''
        node = node.currentText()
        if self.As[self.nodes.index(node)] != 0:
            self.As[self.nodes.index(node)] = eval(origin_A.text())
            state_A.setText("状态: 注意, 该隔室的初始活度已被覆盖")
            return
        self.As[self.nodes.index(node)] = eval(origin_A.text())
        state_A.setText("状态: 添加成功")

    def cal_in(self, time, A_pic):
        '''
        该函数用于计算内照射的剂量
        '''
        A_0 = self.As
        n = len(self.nodes)
        T = self.T
        B = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                B[i][j] = self.edges[i][j] * np.log(2) / self.Ts[j]
                if i == j:
                    B[i][j] = -np.log(2) / self.Ts[j] - np.log(2) / T
        
        def dA_dt(t, A):
            return B @ A
        
        sol = solve_ivp(dA_dt, [0, float(time.text())], A_0, method='RK45', dense_output=True)

        t = np.linspace(0, float(time.text()), 100)
        A = sol.sol(t)

        plt.figure(figsize=(3, 3))
        plt.clf()
        for i in range(n):
            plt.plot(t, A[i], label=self.nodes[i])
        plt.legend()
        plt.savefig("A.png")
        pixmap = QPixmap("A.png")
        pic = QLabel(self)
        pic.setPixmap(pixmap)
        A_pic.setWidget(pic)
