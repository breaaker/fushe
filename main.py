#!/usr/bin/python
import sys
from PyQt6.QtWidgets import QApplication, QMainWindow, QStackedWidget
import cal
import search_formula as sef
import search_params as sep
import ai

class mainwindow(QMainWindow):
    def __init__(self):
        super().__init__()
        
        self.stack = QStackedWidget()
        self.initUI()

    def initUI(self):
        ############################
        #各个页面的初始化
        ############################
        
        self.setWindowTitle("辐射防护应用")

        page_0 = sef.Page_search()
        page_1 = sep.Page_search_gamma()
        page_2 = sep.Page_search_beta()
        page_3 = sep.Page_search_neutron()
        page_4 = cal.Page_cal_out_beta()
        page_5 = cal.Page_cal_out_gamma()
        page_6 = cal.Page_cal_out_neutron()
        page_7 = cal.Page_cal_in()
        page_8 = ai.Page_ai()

        self.stack.addWidget(page_0)
        self.stack.addWidget(page_1)
        self.stack.addWidget(page_2)
        self.stack.addWidget(page_3)
        self.stack.addWidget(page_4)
        self.stack.addWidget(page_5)
        self.stack.addWidget(page_6)
        self.stack.addWidget(page_7)
        self.stack.addWidget(page_8)

        self.stack.setCurrentIndex(0)
        self.setCentralWidget(self.stack)

        ############################
        #定义菜单栏
        ############################

        menu = self.menuBar()
        search_menu = menu.addMenu("查找")
        search_formulas = search_menu.addAction("查公式")
        search_gamma = search_menu.addAction("Gamma相关参数查找")
        search_beta = search_menu.addAction("Beta相关参数查找")
        search_neutron = search_menu.addAction("Neutron相关参数查找")

        cal_menu = menu.addMenu("计算")
        cal_out_beta = cal_menu.addAction("Beta的外照射计算")
        cal_out_gamma = cal_menu.addAction("Gamma的外照射计算")
        cal_out_neutron = cal_menu.addAction("Neutron的外照射计算")
        cal_in = cal_menu.addAction("内照射计算")
        ai_menu = menu.addMenu("辐射防护小助手")
        ai_action = ai_menu.addAction("辐射防护小助手")

        search_formulas.triggered.connect(lambda: self.stack.setCurrentIndex(0))
        search_gamma.triggered.connect(lambda: self.stack.setCurrentIndex(1))
        search_beta.triggered.connect(lambda: self.stack.setCurrentIndex(2))
        search_neutron.triggered.connect(lambda: self.stack.setCurrentIndex(3))
        cal_out_beta.triggered.connect(lambda: self.stack.setCurrentIndex(4))
        cal_out_gamma.triggered.connect(lambda: self.stack.setCurrentIndex(5))
        cal_out_neutron.triggered.connect(lambda: self.stack.setCurrentIndex(6))
        cal_in.triggered.connect(lambda: self.stack.setCurrentIndex(7))
        ai_action.triggered.connect(lambda: self.stack.setCurrentIndex(8))

def main():
    app = QApplication(sys.argv)
    mw = mainwindow()
    mw.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()