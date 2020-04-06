#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPaintEvent>
#include <QPainter>
#include "qcustomplot.h"
#include "Simulation/parameters.hpp"
#include "Simulation/simulation.hpp"
#include "Simulation/random_thijs.hpp"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE
enum display_color {cells, normal_rate, cancer_rate, infected_rate, resistant_rate, dominant_rate};

enum grid_type {regular, voronoi};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
     ~MainWindow();

    void update_image(size_t sq_size);

    void update_image(size_t sq_size,
                      const std::array< binned_distribution, 4 > & growth_rate);


    void display_voronoi(size_t sq_size); // cell coloring
    void display_voronoi(const binned_distribution& growth_rate,
                         cell_type focal_cell_type,
                         size_t sq_size); // growth rate coloring
    void display_voronoi(const std::array< binned_distribution, 4 > & growth_rate,
                         size_t sq_size); // dominant growth rate coloring

    void display_regular(); // cell type coloring
    void display_regular(const binned_distribution& growth_rate,
                         cell_type focal_cell_type); // growth rate coloring
    void display_regular(const std::array< binned_distribution, 4 > & growth_rate); // dominant growth rate coloring


    void update_parameters(Param& p);
    void print_params(const Param& p);

    void set_resolution(int width, int height);
    void set_pixel(int x, int y, const QColor& col);
    void update_polygons(const std::vector< std::vector< voronoi_point > >& all_edges);

    void update_plot(double t, const std::array<size_t, 5>& cell_numbers);
    void setup_simulation();

private slots:
    void on_btn_start_clicked();

    void on_btn_stop_clicked();

    void on_speed_slider_actionTriggered(int action);

    void on_drpdwnbox_display_activated(int index);

    void on_btn_setup_clicked();

    void on_btn_add_virus_clicked();

private:
    Ui::MainWindow *ui;
    QImage image_;

    int row_size;
    int col_size;
    float factor_x;
    float factor_y;

    display_color focal_display_type;

    QVector<double> x_t;
    QVector<double> y_n;
    QVector<double> y_c;
    QVector<double> y_i;
    QVector<double> y_r;

    std::vector< QPolygonF > polygons;

    bool is_running;
    bool is_paused;
    int update_speed;

    grid_type grid_type;

    Param all_parameters;
   // simulation Simulation;
    std::unique_ptr<simulation> sim;

    std::vector< QColor > colorz;
};
#endif // MAINWINDOW_H
