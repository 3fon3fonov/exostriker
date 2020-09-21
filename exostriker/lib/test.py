import pyqtgraph as pg

app = pg.mkQApp()

pw = pg.plot([10, 100, 1000, 100000], [0.1, 0.8, 0.2, 0.4])
pw.plotItem.setLogMode(x=True)

arrow = pg.ArrowItem(pos=(2, 0.8), angle=270)
pw.addItem(arrow)

text = pg.TextItem("hey", anchor=(0.5, 2))
text.setPos(2, 0.8)
pw.addItem(text)

app.exec_()
