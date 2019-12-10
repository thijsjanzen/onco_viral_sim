#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPaintEvent>
#include <QPainter>
#include "qcustomplot.h"
#include "Simulation/node.hpp"
#include "Simulation/parameters.hpp"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE
enum display_color {cells, normal_rate, cancer_rate, infected_rate, resistant_rate, dominant_rate};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
     ~MainWindow();

    void update_image(const std::vector< node >& world, size_t sq_size);
    void update_image(const std::vector< node >& world, size_t sq_size,
                                  const std::vector< std::vector< float> > & growth_rate);

    void update_parameters(Param& p);
    void print_params(const Param& p);

    void set_resolution(int width, int height);
    void set_pixel(int x, int y, const QColor& col);

    void update_plot(double t, const std::vector<int>& cell_numbers);

private slots:
    void on_btn_start_clicked();

    void on_btn_stop_clicked();

    void on_speed_slider_actionTriggered(int action);

    void on_drpdwnbox_display_activated(int index);

private:
    Ui::MainWindow *ui;
    QImage image_;
    display_color focal_display_type;

    QVector<double> x_t;
    QVector<double> y_n;
    QVector<double> y_c;
    QVector<double> y_i;
    QVector<double> y_r;

    bool is_running;
    int update_speed;
};
#endif // MAINWINDOW_H
