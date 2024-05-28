from PyQt6.QtWidgets import (QWidget, QPushButton, QHBoxLayout, QVBoxLayout, QLabel, QLineEdit, QComboBox)
import numpy as np

with open("data/gamma_data/point.txt", "r") as f:
    data = f.readlines()
gamma_point = {}
for i in range(len(data)):
    gamma_data = data[i].strip().split("\t")
    gamma_point[gamma_data[0]] = [float(gamma_data[1]), float(gamma_data[2])]

gamma_fm = np.loadtxt("data/gamma_data/fm.txt", delimiter="\t", skiprows=1, encoding="utf-8")
gamma_fx = np.loadtxt("data/gamma_data/fx.txt", delimiter="\t")

class Page_search_gamma(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):

        ############################
        #gamma
        ############################
        
        #照射量率常数和空气比释动能率常数
        element = QLineEdit(self)
        element.setPlaceholderText("请输入元素名称:例如Na-24")
        Gamma = QLabel()
        Gamma.setText("照射量率常数")
        Gamma_k = QLabel()
        Gamma_k.setText("空气比释动能率常数")
        button = QPushButton("Search", self)
        button.clicked.connect(lambda: self.search_gamma(element, Gamma, Gamma_k))

        layout = QHBoxLayout()#水平布局
        self.setLayout(layout)

        #垂直布局
        se_layout = QVBoxLayout()
        se_layout.addWidget(QLabel("#############################################"))
        se_layout.addWidget(QLabel("照射量率常数和空气比释动能率常数"))
        se_layout.addWidget(element)
        se_layout.addWidget(button)
        se_layout.addWidget(Gamma)
        se_layout.addWidget(Gamma_k)

        #照射量到吸收剂量的转换因子
        
        factor = QComboBox(self)
        factor.addItem("骨")
        factor.addItem("肌肉")
        factor.addItem("软组织")
        factor.addItem("水")

        energy_1 = QLineEdit(self)
        energy_1.setPlaceholderText("请输入gamma能量, 单位MeV, 范围为0.01-10MeV")
        f_m = QLabel()
        f_m.setText("吸收剂量转换因子为：")
        button2 = QPushButton("Search", self)
        button2.clicked.connect(lambda: self.search_f_m(energy_1, factor, f_m))

        se_layout.addWidget(QLabel("#############################################"))
        se_layout.addWidget(QLabel("照射量到吸收剂量的转换因子"))
        se_layout.addWidget(factor)
        se_layout.addWidget(energy_1)
        se_layout.addWidget(button2)
        se_layout.addWidget(f_m)

        se_layout_2 = QVBoxLayout()

        #照射量因子
        energy_2 = QLineEdit(self)
        energy_2.setPlaceholderText("请输入gamma能量, 单位MeV, 范围为0.01-10MeV")
        f_x = QLabel()
        f_x.setText("照射量因子为：")
        button3 = QPushButton("Search", self)
        button3.clicked.connect(lambda: self.search_f_x(energy_2, f_x))
        
        se_layout_2.addWidget(QLabel("#############################################"))
        se_layout_2.addWidget(QLabel("照射量因子"))
        se_layout_2.addWidget(energy_2)
        se_layout_2.addWidget(button3)
        se_layout_2.addWidget(f_x)

        #质量衰减系数， 质量能量吸收系数

        material = QComboBox(self)
        material.addItem("干燥空气")
        material.addItem("水")
        material.addItem("密质骨")
        material.addItem("软组织")
        material.addItem("混凝土")
        material.addItem("有机玻璃")
        material.addItem("元素")
        element_2 = QLineEdit(self)
        element_2.setPlaceholderText("如果选择元素，请输入元素名称")
        energy_3 = QLineEdit(self)
        energy_3.setPlaceholderText("请输入gamma能量, 单位MeV, 范围为0.001-20MeV")
        mu = QLabel()
        mu.setText("质量衰减系数为：")
        mu_en = QLabel()
        mu_en.setText("质量能量吸收系数为：")
        button4 = QPushButton("Search", self)
        button4.clicked.connect(lambda: self.search_mu(material.currentText(), element_2.text(), energy_3.text(), mu, mu_en))

        se_layout_2.addWidget(QLabel("#############################################"))
        se_layout_2.addWidget(QLabel("质量衰减系数， 质量能量吸收系数"))
        se_layout_2.addWidget(material)
        se_layout_2.addWidget(element_2)
        se_layout_2.addWidget(energy_3)
        se_layout_2.addWidget(button4)
        se_layout_2.addWidget(mu)
        se_layout_2.addWidget(mu_en)
        layout.addLayout(se_layout)
        layout.addLayout(se_layout_2)

    def search_f_m(self, energy_1, factor, f_m):
        '''
        该函数用于搜索gamma照射量到吸收剂量的转换因子
        '''
        factor = factor.currentText()
        if energy_1.text() == "":
            f_m.setText("请输入能量")
            return
        energy_1 = float(energy_1.text())
        if factor == "骨":
            f_m.setText("吸收剂量转换因子为：" + linear(energy_1, gamma_fm[:, 0], gamma_fm[:, 1]) + "J/C")
        elif factor == "肌肉":
            f_m.setText("吸收剂量转换因子为：" + linear(energy_1, gamma_fm[:, 0], gamma_fm[:, 2]) + "J/C")
        elif factor == "软组织":
            f_m.setText("吸收剂量转换因子为：" + linear(energy_1, gamma_fm[:, 0], gamma_fm[:, 3]) + "J/C")
        elif factor == "水":
            f_m.setText("吸收剂量转换因子为：" + linear(energy_1, gamma_fm[:, 0], gamma_fm[:, 4]) + "J/C")

    def search_gamma(self, element, Gamma, Gamma_k):
        '''
        该函数用于搜索gamma照射量率常数和空气比释动能率常数
        '''
        element = element.text()
        if element not in gamma_point.keys():
            Gamma.setText("未找到元素")
            Gamma_k.setText("未找到元素")
        else:
            Gamma.setText("照射量率常数: " + str(gamma_point[element][0]) + "Cm^2/kg")
            Gamma_k.setText("空气比释动能率常数: " + str(gamma_point[element][1]) + "Gym^2")

    def search_f_x(self, energy_2, f_x):
        '''
        该函数用于搜索gamma照射量因子
        '''
        if energy_2.text() == "":
            f_x.setText("请输入能量")
            return
        energy_2 = float(energy_2.text())
        f_x.setText("照射量因子为：" + linear(energy_2, gamma_fx[:, 0], gamma_fx[:, 1]) + "Cm^2/kg")

    def search_mu(self, material, element_2, energy_3, mu, mu_en):
        if material == "元素":
            if element_2 == "":
                mu.setText("请输入元素名称")
                mu_en.setText("请输入元素名称")
                return
            if energy_3 == "":
                mu.setText("请输入能量")
                mu_en.setText("请输入能量")
                return
            data = np.loadtxt("data/gamma_data/mu_en_rho_data/" + element_2 + ".txt", delimiter=",")
            mu.setText("质量衰减系数为：" + linear(float(energy_3), data[:, 0], data[:, 1]) + "cm^2/g")
            mu_en.setText("质量能量吸收系数为：" + linear(float(energy_3), data[:, 0], data[:, 2]) + "cm^2/g")
        else:
            if energy_3 == "":
                mu.setText("请输入能量")
                mu_en.setText("请输入能量")
                return
            data = np.loadtxt("data/gamma_data/mu_en_rho_data/" + material + ".txt", delimiter="\t")
            mu.setText("质量衰减系数为：" + linear(float(energy_3), data[:, 0]/1e6, data[:, 1]) + "cm^2/g")
            mu_en.setText("质量能量吸收系数为：" + linear(float(energy_3), data[:, 0]/1e6, data[:, 2]) + "cm^2/g")

phi_0 = np.loadtxt("data/beta_data/单能电子.txt", delimiter="\t", skiprows=1, encoding="utf-8")
with open("data/beta_data/source.txt", "r", encoding="utf-8") as f:
    data = f.readlines()
beta_sources = {}
for i in range(len(data)):
    beta_data = data[i].strip().split("\t")
    beta_sources[beta_data[0]] = [beta_data[1], beta_data[2], beta_data[3]]

class Page_search_beta(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        ############################
        #beta
        ############################
        
        #相当于剂量率为2.5e-3mGy/h的单能电子的注量率值
        energy = QLineEdit(self)
        energy.setPlaceholderText("请输入电子能量, 单位MeV, 范围为0.01-10MeV")
        phi = QLabel()
        phi.setText("相当于剂量率为2.5e-3mGy/h的单能电子的注量率值为: ")
        button = QPushButton("Search", self)
        button.clicked.connect(lambda: self.search_phi(energy, phi))

        layout = QHBoxLayout()#水平布局
        self.setLayout(layout)

        se_layout = QVBoxLayout()
        se_layout.addWidget(QLabel("#############################################"))
        se_layout.addWidget(QLabel("单能电子的注量率值"))
        se_layout.addWidget(energy)
        se_layout.addWidget(button)
        se_layout.addWidget(phi)

        #beta源特性
        source = QLineEdit(self)
        source.setPlaceholderText("请输入源名称, 例如Sr-90等核素, 还有混合裂变产物、天然铀")
        source_info = QLabel()
        source_info.setText("beta源特性为: ")
        button2 = QPushButton("Search", self)
        button2.clicked.connect(lambda: self.search_source(source, source_info))

        se_layout.addWidget(QLabel("#############################################"))
        se_layout.addWidget(QLabel("beta源特性"))
        se_layout.addWidget(source)
        se_layout.addWidget(button2)
        se_layout.addWidget(source_info)

        #质量阻止本领和CSDA射程
        material = QComboBox(self)
        material.addItem("干燥空气")
        material.addItem("水")
        material.addItem("密质骨")
        material.addItem("组织等效塑料")
        material.addItem("石墨")
        material.addItem("眼睛晶状体")
        material.addItem("元素")

        element = QLineEdit(self)
        element.setPlaceholderText("如果选择元素，请输入元素名称")
        energy_2 = QLineEdit(self)
        energy_2.setPlaceholderText("请输入电子能量, 单位MeV, 范围为0.01-1000MeV")
        S_rho = QLabel()
        S_rho.setText("质量阻止本领为：")
        beta_range = QLabel()
        beta_range.setText("CSDA射程为: ")
        button3 = QPushButton("Search", self)
        button3.clicked.connect(lambda: self.search_S_rho(material.currentText(), element.text(), energy_2.text(), S_rho, beta_range))

        se_layout_2 = QVBoxLayout()

        se_layout_2.addWidget(QLabel("#############################################"))
        se_layout_2.addWidget(QLabel("质量阻止本领和CSDA射程"))
        se_layout_2.addWidget(material)
        se_layout_2.addWidget(element)
        se_layout_2.addWidget(energy_2)
        se_layout_2.addWidget(button3)
        se_layout_2.addWidget(S_rho)
        se_layout_2.addWidget(beta_range)

        layout.addLayout(se_layout)
        layout.addLayout(se_layout_2)

    def search_phi(self, energy, phi, cal=False):
        '''
        该函数用于搜索beta相当于剂量率为2.5e-3mGy/h的单能电子的注量率值
        '''
        if energy.text() == "":
            phi.setText("请输入能量")
            return
        energy = float(energy.text())
        if cal:
            return linear(energy, phi_0[:, 0], phi_0[:, 1]*10000)
        phi.setText("相当于剂量率为2.5e-3mGy/h的单能电子的注量率值为: " + linear(energy, phi_0[:, 0], phi_0[:, 1]*10000) + "/(m^2s)")
    
    def search_source(self, source, source_info, cal=False):
        '''
        该函数用于搜索beta源特性
        '''
        source = source.text()
        if source not in beta_sources.keys():
            if cal:
                return [None, None, None]
            source_info.setText("未找到源")
        else:
            if cal:
                return beta_sources[source]
            source_info.setText("beta源特性为: \n半衰期: " + beta_sources[source][0] + ", \nbeta射线最大能量: " + beta_sources[source][1] + ", \nbeta射线平均能量: " + beta_sources[source][2])
    
    def search_S_rho(self, material, element, energy, S_rho, beta_range):
        '''
        该函数用于搜索beta质量阻止本领和CSDA射程
        '''
        if material == "元素":
            if element == "":
                S_rho.setText("请输入元素名称")
                beta_range.setText("请输入元素名称")
                return
            if energy == "":
                S_rho.setText("请输入能量")
                beta_range.setText("请输入能量")
                return
            data = np.loadtxt("data/beta_data/S_rho_range_data/" + element + ".txt")
            S_rho.setText("质量阻止本领为: " + linear(float(energy), data[:, 0], data[:, 1]) + "MeV/cm^2")
            beta_range.setText("CSDA射程为: " + linear(float(energy), data[:, 0], data[:, 2]) + "g/cm^2")
        else:
            if energy == "":
                S_rho.setText("请输入能量")
                beta_range.setText("请输入能量")
                return
            data = np.loadtxt("data/beta_data/S_rho_range_data/" + material + ".txt")
            S_rho.setText("质量阻止本领为: " + linear(float(energy), data[:, 0], data[:, 1]) + "MeV/cm^2")
            beta_range.setText("CSDA射程为: " + linear(float(energy), data[:, 0], data[:, 2]) + "g/cm^2")

with open("data/n_data/sigma.txt", "r", encoding="utf-8") as f:
    data = f.readlines()
n_sigma = {}
for i in range(len(data)):
    n_data = data[i].strip().split("\t")
    n_sigma[n_data[0]] = float(n_data[1])

with open("data/n_data/source.txt", "r", encoding="utf-8") as f:
    data = f.readlines()
n_sources = {}
for i in range(len(data)):
    n_data = data[i].strip().split("\t")
    n_sources[n_data[0]] = [n_data[1], n_data[2], n_data[3], n_data[4], n_data[5]]

class Page_search_neutron(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
    def initUI(self):
        ############################
        #neutron
        ############################

        #比释动能因子f_k
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
        f_k = QLabel()
        f_k.setText("比释动能因子为：")
        button = QPushButton("Search", self)
        button.clicked.connect(lambda: self.search_f_k(material.currentText(), energy.text(), f_k))

        layout = QHBoxLayout()#水平布局
        self.setLayout(layout)

        se_layout = QVBoxLayout()
        se_layout.addWidget(QLabel("Neutron部分"))
        se_layout.addWidget(QLabel("#############################################"))
        se_layout.addWidget(QLabel("比释动能因子"))
        se_layout.addWidget(material)
        se_layout.addWidget(energy)
        se_layout.addWidget(button)
        se_layout.addWidget(f_k)

        #宏观分出截面
        material_2 = QComboBox(self)
        material_2.addItem("普通土(含水10%)")
        material_2.addItem("石墨(密度为1.54)")
        material_2.addItem("普通混凝土")
        material_2.addItem("水")
        material_2.addItem("石蜡")
        material_2.addItem("聚乙烯")
        material_2.addItem("铁")
        material_2.addItem("铅")
        material_2.addItem("铝")
        sigma = QLabel()
        sigma.setText("宏观分出截面为：")
        button2 = QPushButton("Search", self)
        button2.clicked.connect(lambda: self.search_n_sigma(material_2.currentText(), sigma))

        se_layout.addWidget(QLabel("#############################################"))
        se_layout.addWidget(QLabel("宏观分出截面"))
        se_layout.addWidget(material_2)
        se_layout.addWidget(button2)
        se_layout.addWidget(sigma)

        #中子源特性
        n_source = QComboBox(self)
        n_source.addItem("Na-Be")
        n_source.addItem("Sb-Be")
        n_source.addItem("Po-Be")
        n_source.addItem("Ra-Be")
        n_source.addItem("Pu8-Be")
        n_source.addItem("Pu9-Be")
        n_source.addItem("Am-Be")

        source_info = QLabel()
        source_info.setText("中子源特性为: ")
        button3 = QPushButton("Search", self)
        button3.clicked.connect(lambda: self.search_n_source(n_source, source_info))

        se_layout_2 = QVBoxLayout()
        se_layout_2.addWidget(QLabel("#############################################"))
        se_layout_2.addWidget(QLabel("中子源特性"))
        se_layout_2.addWidget(n_source)
        se_layout_2.addWidget(button3)
        se_layout_2.addWidget(source_info)
    
        #品质因子, 剂量当量换算因子, 中子注量率
        energy_2 = QLineEdit(self)
        energy_2.setPlaceholderText("请输入中子能量, 单位MeV, 范围为2.5e-8到50MeV")
        Q = QLabel()
        Q.setText("品质因子为：")
        f_n = QLabel()
        f_n.setText("剂量当量换算因子为：")
        phi = QLabel()
        phi.setText("与10muSvm^2/h相对应的中子注量率为: ")

        button4 = QPushButton("Search", self)
        button4.clicked.connect(lambda: self.search_Q_f_n_phi(energy_2.text(), Q, f_n, phi))

        se_layout_2.addWidget(QLabel("#############################################"))
        se_layout_2.addWidget(energy_2)
        se_layout_2.addWidget(button4)
        se_layout_2.addWidget(Q)
        se_layout_2.addWidget(f_n)
        se_layout_2.addWidget(phi)

        layout.addLayout(se_layout)
        layout.addLayout(se_layout_2)

    def search_f_k(self, material, energy, f_k, cal=False):
        if energy == "":
            f_k.setText("请输入能量")
            return
        energy = float(energy)

        data = np.loadtxt("data/n_data/fk_data/" + material + ".txt")
        if cal:
            return float(linear(energy, data[:, 0], data[:, 1]))
        f_k.setText("比释动能因子为：" + linear(energy, data[:, 0], data[:, 1]) + " Gycm^2")

    def search_n_sigma(self, material, sigma, cal=False):
        if cal:
            return n_sigma[material]
        sigma.setText("宏观分出截面为：" + str(n_sigma[material]) + "/cm")

    def search_n_source(self, n_source, source_info, cal=False):
        source = n_source.currentText()
        if source not in n_sources.keys():
            if cal:
                return [None, None, None, None]
            source_info.setText("未找到源")
        else:
            if cal:
                return n_sources[source][1:]
            source_info.setText("中子源特性为: \n半衰期: " + n_sources[source][0] + ", \n中子最大能量: " + n_sources[source][1] + "MeV, \n中子平均能量: " + n_sources[source][2] + "MeV, \n中子产额: " + n_sources[source][3])
    
    def search_Q_f_n_phi(self, energy_2, Q, f_n, phi):
        if energy_2 == "":
            Q.setText("请输入能量")
            f_n.setText("请输入能量")
            phi.setText("请输入能量")
            return
        energy_2 = float(energy_2)
        data = np.loadtxt("data/n_data/Q_fn_phi.txt")
        Q.setText("品质因子为：" + linear(energy_2, data[:, 0], data[:, 1]))
        f_n.setText("剂量当量换算因子为：" + linear(energy_2, data[:, 0], data[:, 2]) + "e-15Svm^2")
        phi.setText("与10muSvm^2/h相对应的中子注量率为: " + linear(energy_2, data[:, 0], data[:, 3]) + "/(cm^2s)")

def linear(energy, v_1, v_2):
    if energy < v_1[0]:
        return "energy too low"
    elif energy > v_1[-1]:
        return "energy too high"
    else:
        return str(np.interp(energy, v_1, v_2))