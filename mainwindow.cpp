#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QPaintEvent>
#include <QPainter>
#include <string>
#include <sstream>
#include <chrono>
#include <iostream>

#include "Simulation/parameters.hpp"
#include "Simulation/simulation.hpp"
#include "Simulation/rndutils.hpp"


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->line_plot->addGraph(); // normal
    ui->line_plot->addGraph(); // cancer
    ui->line_plot->addGraph(); // infected
    ui->line_plot->graph(0)->setPen(QPen(Qt::blue));
    ui->line_plot->graph(1)->setPen(QPen(Qt::red));
    ui->line_plot->graph(2)->setPen(QPen(Qt::green));

    ui->line_plot->graph(0)->setName("Normal");
    ui->line_plot->graph(1)->setName("Cancer");
    ui->line_plot->graph(2)->setName("Infected");

    QCPPlotTitle *fst_title = new QCPPlotTitle(ui->line_plot, "Number of cells");
    ui->line_plot->plotLayout()->insertRow(0);
    ui->line_plot->plotLayout()->addElement(0, 0, fst_title);
    ui->line_plot->xAxis->setLabel("Time (hours)");
    ui->line_plot->yAxis->setLabel("Number of Cells");

    ui->line_plot->legend->setVisible(true);
    QFont legendFont = font();
    legendFont.setPointSize(7);
    ui->line_plot->legend->setFont(legendFont);
    ui->line_plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignLeft);

    ui->progressBar->setValue(0);

    update_speed = ui->speed_slider->value();

    ui->box_infection_routine->addItem("Center");
    ui->box_infection_routine->addItem("Random");

    ui->box_start_setup->addItem("Full");
    ui->box_start_setup->addItem("Grow");
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::set_resolution(int width, int height) {
    image_ = QImage(width,height, QImage::Format_RGB32);
}

void MainWindow::set_pixel(int x, int y, const QColor& col) {
    image_.setPixel(x, y, col.rgb());
}


void MainWindow::update_image(const std::vector< node >& world) {
    QPainter p(&image_);

    std::vector< QColor > colorz = {
        {0, 0, 255}, // blue, normal
        {255, 0, 0}, // red,  cancer
        {0, 255, 0}, // green, infected
        {0, 0, 0} // black, empty
    };

    for(auto i : world) {
      QColor c = colorz[i.node_type];
      p.setPen(c);
      p.drawPoint(static_cast<int>(i.x_),
                  static_cast<int>(i.y_));
    }

    int w = ui->q_label->width();
    int h = ui->q_label->height();

     ui->q_label->setPixmap((QPixmap::fromImage(image_)).scaled(w,h, Qt::KeepAspectRatio));
     ui->q_label->repaint();
     ui->q_label->update();
     ui->q_label->show();
     QApplication::processEvents();
}

std::string get_string(std::string s, float v) {
    std::string output = s + " " + std::to_string(v) + "\n";
    return output;
}

void MainWindow::print_params(const Param& p) {
    // TODO
    std::stringstream s;
    s << get_string("Maximum Time", p.maximum_time);
    s << get_string("Cancer time", p.time_adding_cancer);
    s << get_string("Virus time", p.time_adding_cancer);
    s << get_string("Initial # cancer cells", p.initial_number_cancer_cells);
    s << get_string("Initial # of normal cells", p.initial_number_normal_cells);
    s << get_string("Birth Rate Normal", p.birth_normal);
    s << get_string("Death Rate Normal", p.death_normal);
    s << get_string("Birth Rate Cancer", p.birth_cancer);
    s << get_string("Death Rate Cancer", p.death_cancer);
    s << get_string("Birth Rate Infected", p.birth_infected);
    s << get_string("Death Rate Infected", p.death_infected);
    s << get_string("Infection routine", p.infection_type);
    s << get_string("Infection %", p.percent_infected);
    s << get_string("Start type", p.start_setup);

    ui->text->appendPlainText(QString::fromStdString(s.str()));
    return;
}


void MainWindow::update_parameters(Param& p) {

   p.maximum_time = static_cast<int>(ui->box_maxtime->value());
   p.time_adding_cancer = static_cast<int>(ui->box_cancer_time->value());
   p.time_adding_virus = static_cast<int>(ui->box_virus_time->value());

   p.initial_number_cancer_cells = static_cast<int>(ui->box_cancer_cells->value());
   p.initial_number_normal_cells = static_cast<int>(ui->box_normal_cells->value());

   p.birth_normal = static_cast<float>(ui->box_birth_normal->value());
   p.death_normal = static_cast<float>(ui->box_death_normal->value());

   p.birth_cancer = static_cast<float>(ui->box_birth_cancer->value());
   p.death_cancer = static_cast<float>(ui->box_death_cancer->value());

   p.birth_infected = static_cast<float>(ui->box_birth_infected->value());
   p.death_infected = static_cast<float>(ui->box_death_infected->value());

   p.percent_infected = static_cast<float>(ui->box_percent_infected->value());
   p.infection_type = random_infection;

   auto infection_string = ui->box_infection_routine->currentText();
   if(infection_string == "Random")
       p.infection_type = random_infection;
   if(infection_string == "Center")
       p.infection_type = center_infection;


   auto start_string = ui->box_start_setup->currentText();
    if(start_string == "Grow")
        p.start_setup = grow;
    if(start_string == "Full")
        p.start_setup = full;

   print_params(p);
   return;
}


void MainWindow::on_btn_start_clicked()
{
    auto start_time = std::chrono::high_resolution_clock::now();

    x_t.clear();
    y_n.clear();
    y_c.clear();
    y_i.clear();

    Param all_parameters;
    update_parameters(all_parameters);

    simulation Simulation(all_parameters);

    Simulation.initialize_network();

    set_resolution(static_cast<int>(Simulation.sq_size),
                   static_cast<int>(Simulation.sq_size));

    Simulation.t = 0.0;
    int counter = 0;
    is_running = true;



    while(Simulation.t < all_parameters.maximum_time) {
        Simulation.update_one_step();
        counter++;

        int progress = static_cast<int>(100.f * Simulation.t / all_parameters.maximum_time);
        ui->progressBar->setValue(progress);
        int update_step = static_cast<int>((1.0f * update_speed / 100) * Simulation.world.size());
        if(counter % update_step == 0) {
            update_image(Simulation.world);

            update_plot(static_cast<double>(Simulation.t),
                        Simulation.get_cell_numbers());
        }
        if(!is_running) break;
    }


    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( end_time - start_time ).count();

    float time_taken = 1.f * duration / 1000;

    std::stringstream st;
    st << "This took ";
    st << time_taken;
    st << " seconds\n";
    ui->text->appendPlainText(QString::fromStdString(st.str()));
    std::cout << time_taken<< "\n";

}

void MainWindow::update_plot(double t, const std::vector<int> &cell_numbers) {
    x_t.append(t);
    y_n.append(cell_numbers[0]);
    y_c.append(cell_numbers[1]);
    y_i.append(cell_numbers[2]);

    ui->line_plot->graph(0)->clearData();
    ui->line_plot->graph(0)->setData(x_t, y_n);

    ui->line_plot->graph(1)->clearData();
    ui->line_plot->graph(1)->setData(x_t, y_c);

    ui->line_plot->graph(2)->clearData();
    ui->line_plot->graph(2)->setData(x_t, y_i);

    ui->line_plot->rescaleAxes();
    ui->line_plot->replot();
    ui->line_plot->update();
    QApplication::processEvents();
}



void MainWindow::on_btn_stop_clicked() {
    is_running = false;
}

void MainWindow::on_speed_slider_actionTriggered(int action) {
   update_speed = ui->speed_slider->value();
}
