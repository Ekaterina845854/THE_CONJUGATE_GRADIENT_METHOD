#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>
#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>
#include <QDebug>
#include <QVBoxLayout>
#include <QFile>
#include <QTextStream>
#include <QDebug>
#include <QDir>

using namespace std;
using namespace QtCharts;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow) {
    ui->setupUi(this);

    if (!ui->chartView) {
        qDebug() << "Error: chartView is not initialized!";
        return;
    }
}

MainWindow::~MainWindow() {
    delete ui;
}

// Определение метода vectorAdd
std::vector<double> MainWindow::vectorAdd(const std::vector<double>& v1, const std::vector<double>& v2) {
    std::vector<double> result;
    for (size_t i = 0; i < v1.size(); ++i) {
        result.push_back(v1[i] + v2[i]);
    }
    return result;
}

// Определение метода F1
double MainWindow::F1(const std::vector<double>& x) {
    return ((x[0] - 5) * (x[0] - 5) * (x[1] - 4) * (x[1] - 4)) + (20 * (x[0] - 5) * (x[0] - 5)) + ((x[1] - 10) * (x[1] - 10));
}

// Определение метода F2
double MainWindow::F2(const std::vector<double>& x) {
    return -1.0 / (((x[0] - 2) * (x[0] - 2) * (x[1] - 2) * (x[1] - 2)) + (40 * (x[0] - 5) * (x[0] - 5)) + ((x[1] - 4) * (x[1] - 4)) + 1);
}

// Определение метода GradF1
std::vector<double> MainWindow::GradF1(std::function<double(const std::vector<double>&)> f, const std::vector<double>& x) {
    const double h = 1e-5;
    std::vector<double> grad(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        std::vector<double> xh = x;
        xh[i] += h;
        grad[i] = (f(xh) - f(x)) / h;
    }
    return grad;
}

// Определение метода FindMin
double MainWindow::FindMin(std::function<double(const std::vector<double>&)> f, const std::vector<double>& start, const std::vector<double>& direction) {
    double alpha = 0.0;
    double step = 0.1;
    std::vector<double> x = start;
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] += alpha * direction[i];
    }
    double fx = f(x);
    while (true) {
        alpha += step;
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = start[i] + alpha * direction[i];
        }
        double new_fx = f(x);
        if (new_fx < fx) {
            fx = new_fx;
        } else {
            break;
        }
    }
    return fx;
}

// Функция для сохранения результатов в файл
void MainWindow::saveResultsToFile(const QString& fileName, const std::vector<double>& result, int iterations) {
    QString path = "D:/CursachProga/Try_1/" + fileName; // Путь к файлу
    QFile file(path);
    if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream out(&file);
        out << "Result: (" << result[0] << ", " << result[1] << ")\n";
        out << "Iterations: " << iterations << "\n";
        file.close();
        qDebug() << "Результаты успешно сохранены в файл:" << path;
    } else {
        qDebug() << "Ошибка при открытии файла:" << path;
    }
}

// Функция создания файла
void MainWindow::createFile(const QString &fileName) {
    QString directory = "D:\\CursachProga\\Try_1";
    QDir dir(directory);

    // Создаем папку, если она не существует
    if (!dir.exists()) {
        dir.mkpath(directory);
    }

    QFile file(dir.filePath(fileName));
    if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream out(&file);
        out << "Файл создан успешно!\n";  // Можно записать любую информацию
        file.close();
        qDebug() << "Файл создан: " << dir.filePath(fileName);
    } else {
        qDebug() << "Ошибка при создании файла: " << file.errorString();
    }
}


std::vector<double> MainWindow::FletcherReevesMethod(const QString& fileName, std::function<double(const std::vector<double>&)> func, double tolerance, double lambdaTolerance, std::vector<double> initialPoint) {
    std::vector<double> x = initialPoint;
    std::vector<double> p = GradF1(func, x);
    double gradSquare = inner_product(p.begin(), p.end(), p.begin(), 0.0);
    int iteration = 0;

    QLineSeries *series = new QLineSeries();
    series->append(x[0], x[1]);

    double min_x = x[0], max_x = x[0];
    double min_y = x[1], max_y = x[1];

    while (gradSquare > tolerance * tolerance) {
        iteration++;
        double alpha = FindMin(func, x, p);

        for (size_t i = 0; i < x.size(); ++i) {
            x[i] += alpha * p[i];
        }

        // Обновление минимальных и максимальных значений
        min_x = std::min(min_x, x[0]);
        max_x = std::max(max_x, x[0]);
        min_y = std::min(min_y, x[1]);
        max_y = std::max(max_y, x[1]);

        std::vector<double> gradNew = GradF1(func, x);
        double gradNewSquare = inner_product(gradNew.begin(), gradNew.end(), gradNew.begin(), 0.0);
        double beta = gradNewSquare / gradSquare;

        for (size_t i = 0; i < p.size(); ++i) {
            p[i] = gradNew[i] + beta * p[i];
        }

        gradSquare = gradNewSquare;
        series->append(x[0], x[1]);
        qDebug() << "Итерация " << iteration << ": (" << x[0] << ", " << x[1] << ")";
    }

    qDebug() << "Количество итераций: " << iteration;
    visualizeGraph(series, min_x, max_x, min_y, max_y);

    saveResultsToFile(fileName, x, iteration);
    return x;
}

void MainWindow::visualizeGraph(QtCharts::QLineSeries *series, double min_x, double max_x, double min_y, double max_y) {
    QtCharts::QChart *chart = new QtCharts::QChart();
    chart->addSeries(series);
    chart->setTitle("График");

    // Создание осей
    auto xAxis = new QtCharts::QValueAxis();
    auto yAxis = new QtCharts::QValueAxis();

    // Установка диапазонов осей
    xAxis->setRange(min_x - 0.05 * (max_x - min_x), max_x + 0.05 * (max_x - min_x));
    yAxis->setRange(min_y - 0.05 * (max_y - min_y), max_y + 0.05 * (max_y - min_y));

    // Установка формата меток
    xAxis->setLabelFormat("%g");
    yAxis->setLabelFormat("%7.2g");

    // Установка количества делений
    xAxis->setTickCount(5);
    yAxis->setTickCount(5);

    // Добавление осей в график
    chart->setAxisX(xAxis, series);
    chart->setAxisY(yAxis, series);

    // Установка названий осей
    xAxis->setTitleText("Ось X");  // Название оси X
    yAxis->setTitleText("Ось Y");  // Название оси Y

    // Настройка отображения графика
    QtCharts::QChartView *chartView = new QtCharts::QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);

    // Добавление графика в layout
    QVBoxLayout *layout = dynamic_cast<QVBoxLayout*>(ui->chartView->layout());
    if (!layout) {
        layout = new QVBoxLayout(ui->chartView);
        ui->chartView->setLayout(layout);
    }
    layout->addWidget(chartView);
    chart->legend()->setVisible(false);
}

void MainWindow::on_startButton_clicked() {
    double tolerance = ui->toleranceInput->text().toDouble();
    double lambdaTolerance = ui->lambdaToleranceInput->text().toDouble();
    double funcNumber = ui->funcNumber->text().toDouble(); // Получаем номер функции
    QString fileName = ui->fileNameInput->text().isEmpty() ? "output.txt" : ui->fileNameInput->text();

    if (ui->x0Input->text().isEmpty() || ui->y0Input_2->text().isEmpty()) {
        qDebug() << "Ошибка: начальные координаты не заданы.";
        return;
    }

    vector<double> initialPoint = { ui->x0Input->text().toDouble(), ui->y0Input_2->text().toDouble() };
    qDebug() << "Начальные координаты: (" << initialPoint[0] << ", " << initialPoint[1] << ")";
    qDebug() << "Точность: " << tolerance << ", Точность для alpha: " << lambdaTolerance;

    createFile(fileName);

    try {
        // Выбор функции в зависимости от значения funcNumber
        std::function<double(const std::vector<double>&)> selectedFunction;
        if (funcNumber == 1) {
            selectedFunction = std::bind(&MainWindow::F1, this, std::placeholders::_1);
        } else if (funcNumber == 2) {
            selectedFunction = std::bind(&MainWindow::F2, this, std::placeholders::_1);
        } else {
            qDebug() << "Ошибка: Неверный номер функции. Доступные номера: 1 или 2.";
            return;
        }

        vector<double> result = FletcherReevesMethod(fileName, selectedFunction, tolerance, lambdaTolerance, initialPoint);
        qDebug() << "Результат: (" << result[0] << ", " << result[1] << ")";

    } catch (const std::exception& e) {
        qDebug() << "Ошибка: " << e.what();
    }

}


