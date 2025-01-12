// MainWindow.h - Ensure correct declaration of member functions
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtCharts/QLineSeries>
#include <QtCharts/QChartView>
#include <QVBoxLayout>
#include <vector>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

using namespace std;

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    void createFile(const QString &fileName);
    void saveResultsToFile(const QString& fileName, const std::vector<double>& result, int iterations);
    std::vector<double> vectorAdd(const std::vector<double>& v1, const std::vector<double>& v2);
    double F1(const std::vector<double>& x);
    double F2(const std::vector<double>& x) ;
    std::vector<double> GradF1(std::function<double(const std::vector<double>&)> f, const std::vector<double>& x);
    double MainWindow::FindMin(std::function<double(const std::vector<double>&)> f, const std::vector<double>& start, const std::vector<double>& direction);
    ~MainWindow();

private slots:
    void on_startButton_clicked();

private:
    Ui::MainWindow *ui;
    void MainWindow::visualizeGraph(QtCharts::QLineSeries *series, double min_x, double max_x, double min_y, double max_y);
    std::vector<double> FletcherReevesMethod(const QString& fileName, std::function<double(const std::vector<double>&)> func, double tolerance, double lambdaTolerance, std::vector<double> initialPoint);
};



#endif // MAINWINDOW_H
