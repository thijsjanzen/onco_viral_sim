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
enum display_color {cells, inflammation_factor, added_death_rate,
                    normal_rate, normal_death_rate,
                    cancer_rate, cancer_death_rate,
                    infected_rate, resistant_rate,
                    dominant_rate};



class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
     ~MainWindow();

    void update_image(size_t sq_size,
                      int display_t_cells);

    void update_image(size_t sq_size,
                      const std::array< binned_distribution, 4 > & growth_rate);


    void display_voronoi(size_t sq_size,
                         int display_t_cells); // cell coloring
    void display_voronoi(const binned_distribution& growth_rate,
                         cell_type focal_cell_type,
                         size_t sq_size); // growth rate coloring
    void display_voronoi(const std::array< binned_distribution, 4 > & growth_rate,
                         size_t sq_size); // dominant growth rate coloring

    void display_regular(int display_t_cells); // cell type coloring
    void display_regular(const binned_distribution& growth_rate,
                         cell_type focal_cell_type); // growth rate coloring
    void display_regular(const std::array< binned_distribution, 4 > & growth_rate); // dominant growth rate coloring
    void display_regular_death_rate(const binned_distribution& death_rate,
                            cell_type focal_cell_type,
                            size_t sq_size);


    void update_parameters(Param& p);
    void print_params(const Param& p);

    void set_resolution(int width, int height);
    void set_pixel(int x, int y, const QColor& col);
    void update_polygons(const std::vector< std::vector< voronoi_point > >& all_edges);

    void update_plot(double t, const std::array<size_t, 5>& cell_numbers);
    void setup_simulation();
    void obtain_equilibrium();
    void update_display();

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

    std::unique_ptr<simulation> sim;

    std::vector< QColor > colorz;
};
#endif // MAINWINDOW_H
