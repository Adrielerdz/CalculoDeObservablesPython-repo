import sys 
from PyQt5.QtWidgets import QApplication,QMainWindow 
from PyQt5.QtCore import  QCoreApplication
from lib.design import Ui_Dialog  
import lib.cilindric as cilindric
import lib.gamma as gamma
import lib.tov as tov

class MetricSelectionError(Exception):
    pass

class UnitsSelectionError(Exception):
    pass

class MagnitudIndexError(Exception):
    pass


class MainWindow(Ui_Dialog):

    def setupUi(self,Dialog): 
        super().setupUi(Dialog)
        self.MyActions()

    def MyActions(self):
        self.StartButton.clicked.connect(self.start)
        self.CloseButton.clicked.connect(self.close)

    def close(self):
        QCoreApplication.instance().quit()

    def start(self):
        self.ReadLinesLabel.setText("")
        self.DescartedLinesLabel.setText("")
        self.TimeLabel.setText("")
        
        name_in = self.EoSlineEdit.text()
        name_out = self.OutputFilelineEdit.text()

        try:
            index_E = int(self.EnergyDensitylineEdit.text())
            index_Pl = int(self.ParallelPressureLineEdit.text())
            index_P = int(self.PerpendiculaPressureLineEdit.text())
            index_Rho = int(self.BaryonDensityLineEdit.text())
            index_Nb = int(self.BaryonConcentrationLineEdit.text())
        except ValueError :
            raise MagnitudIndexError

        if self.MevRradioButton.isChecked() and self.MeV3RadioButton.isChecked():
            convertion_factor = 1
        elif self.MeVfm3RadioButton.isChecked() and self.fm3RadioButton.isChecked():
            convertion_factor = 2
        else:
            raise UnitsSelectionError  


        if self.TOVradioButton.isChecked():
            read, descarted , time = tov.main(name_in,name_out,index_E,index_Pl,index_P,index_Rho,index_Nb,convertion_factor)
        elif self.GammaRadioButton.isChecked():
            read, descarted , time = gamma.main(name_in,name_out,index_E,index_Pl,index_P,index_Rho,index_Nb,convertion_factor)
        elif self.CilindricRadioButton.isChecked():
            read, descarted , time = cilindric.main(name_in,name_out,index_E,index_Pl,index_P,index_Rho,index_Nb,convertion_factor)
        else:
            raise  MetricSelectionError


        self.ReadLinesLabel.setText(str(read))
        self.DescartedLinesLabel.setText(str(descarted))
        self.TimeLabel.setText("%d min, %d s"%(int(time)//60, int(time-time//60*60)  ))



if __name__ == '__main__':     
    app = QApplication(sys.argv)     
    w = QMainWindow()    
    ui = MainWindow()     
    ui.setupUi(w)     
    w.show()     
    sys.exit(app.exec_())
