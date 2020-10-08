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

    std::stringstream s;
    s << "Welcome to this In Silico Simulation of oncolytic tumor virotherapy\n";
    s << "Copyright 2019 - 2020 D. Bhatt, T. Janzen & F.J. Weissing\n";
    s << "This is version: 0.5\n";
    ui->text->appendPlainText(QString::fromStdString(s.str()));

    // put the simulation as a unique_ptr
    update_parameters(all_parameters);

    ui->line_plot->addGraph(); // normal
    ui->line_plot->addGraph(); // cancer
    ui->line_plot->addGraph(); // infected
    ui->line_plot->addGraph(); // resistant
    ui->line_plot->graph(0)->setPen(QPen(Qt::blue));
    ui->line_plot->graph(1)->setPen(QPen(Qt::red));
    ui->line_plot->graph(2)->setPen(QPen(Qt::green));
    ui->line_plot->graph(3)->setPen(QPen(Qt::magenta));

    ui->line_plot->graph(0)->setName("Normal");
    ui->line_plot->graph(1)->setName("Cancer");
    ui->line_plot->graph(2)->setName("Infected");
    ui->line_plot->graph(3)->setName("Resistant");

    QCPPlotTitle *fst_title = new QCPPlotTitle(ui->line_plot, "Number of cells");
    ui->line_plot->plotLayout()->insertRow(0);
    ui->line_plot->plotLayout()->addElement(0, 0, fst_title);
    ui->line_plot->xAxis->setLabel("Time (hours)");
    ui->line_plot->yAxis->setLabel("Number of Cells");

    ui->line_plot->legend->setVisible(true);
    QFont legendFont = font();
   // legendFont.setPointSize(7);
    ui->line_plot->legend->setFont(legendFont);
    ui->line_plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignLeft);

    ui->progressBar->setValue(0);

    update_speed = ui->speed_slider->value();

    ui->box_infection_routine->addItem("Center");
    ui->box_infection_routine->addItem("Random");
    ui->box_infection_routine->addItem("Periphery");

    ui->box_start_setup->addItem("Full");
    ui->box_start_setup->addItem("Grow");


    ui->drpdwnbox_display->addItem("Cell types");
    ui->drpdwnbox_display->addItem("T-cells");
    ui->drpdwnbox_display->addItem("Normal Growth Rate");
    ui->drpdwnbox_display->addItem("Normal Death Rate");
    ui->drpdwnbox_display->addItem("Cancer Growth Rate");
    ui->drpdwnbox_display->addItem("Cancer Death Rate");
    ui->drpdwnbox_display->addItem("Infected Growth Rate");
    ui->drpdwnbox_display->addItem("Resistant Growth Rate");
    ui->drpdwnbox_display->addItem("Dominant Growth Rate");

    ui->box_grid_type->addItem("regular");
    ui->box_grid_type->addItem("voronoi");

    is_paused = false;
    is_running = false;

    colorz.push_back(QColor(0, 0, 255));
    colorz.push_back(QColor(255, 0, 0));
    colorz.push_back(QColor(0, 255, 0));
    colorz.push_back(QColor(128, 0, 128));
    colorz.push_back(QColor(0, 0, 0));

    factor_x = 1.f;
    factor_y = 1.f;
}


MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::set_resolution(int width, int height) {
    row_size = width;
    col_size = height;

    image_ = QImage(row_size, col_size, QImage::Format_RGB32);
}

void MainWindow::set_pixel(int x, int y, const QColor& col) {
    image_.setPixel(x, y, col.rgb());
}

QColor get_t_cell_color(float concentration) {
  if(concentration < 1e-5f) {
      QColor col = {0, 0, 0, 255};
      return col;
  }

  QColor col = {255, 0, 255, static_cast<int>(concentration * 255)};
  return col;
}

void MainWindow::display_regular(bool display_t_cells) {
  size_t line_size = row_size;
  size_t num_lines = col_size;

  for(size_t i = 0; i < num_lines; ++i) {
      QRgb* row = (QRgb*) image_.scanLine(i);

      size_t start = i * line_size;
      size_t end = start + line_size;

      for(size_t index = start; index < end; ++index) {
          size_t local_index = index - start;
          if(!display_t_cells) row[local_index] = colorz[ sim->world[index].get_cell_type() ].rgb();
          if(display_t_cells) row[local_index] = get_t_cell_color(sim->world[index].t_cell_concentration).rgba();
      }
  }
}

void MainWindow::display_voronoi(size_t sq_size,
                                 bool display_t_cells) {
  image_.fill(Qt::gray);

  QPainter painter(&image_);

  for(size_t i = 0; i < sim->num_cells; ++i) {
      QBrush brush;
      if(display_t_cells) {
          QColor t_col = get_t_cell_color(sim->world[i].t_cell_concentration);
       /*   float conc = sim->world[i].t_cell_concentration;
          if(conc > 0.01f) {
            t_col.setAlpha(
                static_cast<int>(sim->world[i].t_cell_concentration * 255));
          }*/
          brush = QBrush(t_col);
      } else {
          brush = QBrush(colorz[sim->world[i].get_cell_type()]); // Qt::SolidPattern by default.
      }

      QPainterPath path;
      path.addPolygon(polygons[i]);

      painter.fillPath(path, brush);
      painter.drawPolygon(polygons[i]);
  }
}

void MainWindow::update_image(size_t sq_size,
                              bool display_t_cells) {
    if(grid_type == regular) {
        display_regular(display_t_cells);
    }
    if(grid_type == voronoi) {
        display_voronoi(sq_size, display_t_cells);
    }

    int w = ui->q_label->width();
    int h = ui->q_label->height();

    ui->q_label->setPixmap((QPixmap::fromImage(image_)).scaled(w,h, Qt::KeepAspectRatio));
    ui->q_label->update();
}



QRgb get_color(const cell_type focal_cell_type, float rate) {

    if(rate < 1e-6f) {
        QColor col = {0, 0, 0, 255};
        return col.rgba();
    }

    if(focal_cell_type == normal) {
        QColor col = {0, 0, 255, static_cast<int>(rate * 255)};
        return col.rgba();
    }
    if(focal_cell_type == cancer) {
        QColor col = {255, 0, 0, static_cast<int>(rate * 255)};
        return col.rgba();
    }
    if(focal_cell_type == infected) {
        QColor col = {0, 255, 0, static_cast<int>(rate * 255)};
        return col.rgba();
    }
    if(focal_cell_type == resistant) {
        QColor col = {128, 0, 128, static_cast<int>(rate * 255)};
        return col.rgba();
    }
}

size_t which_max(const std::vector<float>& v) {
    size_t max_index = 0;
    float max_val = v[max_index];
    for(size_t i = 1; i < v.size(); ++i) {
        if(v[i] > max_val) {
            max_index = i;
            max_val = v[i];
        }
    }
    return max_index;
}


void MainWindow::display_regular(const binned_distribution& growth_rate,
                                 cell_type focal_cell_type) {
  int line_size = row_size;
  int num_lines = col_size;

  for(size_t i = 0; i < num_lines; ++i) {
      QRgb* row = (QRgb*) image_.scanLine(i);

      size_t start = i * line_size;
      size_t end = start + line_size;

      for(size_t index = start; index < end; ++index) {
          size_t local_index = index - start;

          row[local_index] = get_color(focal_cell_type, growth_rate.get_value(index));
      }
  }
}

void MainWindow::display_regular(const std::array< binned_distribution, 4 > & growth_rate) {
  size_t line_size = static_cast<size_t>(row_size);
  size_t num_lines = static_cast<size_t>(col_size);

  for(size_t i = 0; i < num_lines; ++i) {
      QRgb* row = (QRgb*) image_.scanLine(i);

      size_t start = i * line_size;
      size_t end = start + line_size;

      for(size_t index = start; index < end; ++index) {
          size_t local_index = index - start;

          std::vector<float> probs = {0.0, 0.0, 0.0, 0.0};
          for(size_t i = 0; i < 4; ++i) {
              probs[i] = growth_rate[i].get_value(index);
          }
          cell_type focal_cell_type = static_cast<cell_type>(which_max(probs));

          row[local_index] = get_color(focal_cell_type,
                                       growth_rate[ focal_cell_type ].get_value(index));
      }
  }
}

void MainWindow::display_voronoi(const binned_distribution& growth_rate,
                                 cell_type focal_cell_type,
                                 size_t sq_size) {
  image_.fill(Qt::gray);

  QPainter painter(&image_);

  for(size_t i = 0; i < sim->world.size(); ++i) {

      QBrush brush(get_color(focal_cell_type,
                             growth_rate.get_value(sim->world[i].pos)));
      QPainterPath path;
      path.addPolygon(polygons[i]);

      painter.fillPath(path, brush);
      painter.drawPolygon(polygons[i]);
  }
}

void MainWindow::display_voronoi(const std::array< binned_distribution, 4 > & growth_rate,
                                 size_t sq_size) {
  image_.fill(Qt::gray);

  QPainter painter(&image_);

  for(size_t i = 0; i < sim->world.size(); ++i) {
      std::vector<float> probs = {0.0, 0.0, 0.0, 0.0};
      for(size_t k = 0; k < 4; ++k) {
          probs[k] = growth_rate[k].get_value(i);
      }
      cell_type focal_cell_type = static_cast<cell_type>(which_max(probs));

      QBrush brush(get_color(focal_cell_type,
                             probs[focal_cell_type]));
      QPainterPath path;
      path.addPolygon(polygons[i]);

      painter.fillPath(path, brush);
      painter.drawPolygon(polygons[i]);
  }
}

void MainWindow::display_regular_death_rate(const binned_distribution& death_rate,
                                    cell_type focal_cell_type,
                                    size_t sq_size) {
  size_t line_size = row_size;
  size_t num_lines = col_size;

  float expected_max = sim->calc_max_t_cell_rate();

  for(size_t i = 0; i < num_lines; ++i) {
      QRgb* row = (QRgb*) image_.scanLine(i);

      size_t start = i * line_size;
      size_t end = start + line_size;

      for(size_t index = start; index < end; ++index) {
          size_t local_index = index - start;

          if(focal_cell_type != cancer) {
              row[local_index] = get_color(focal_cell_type,
                                       death_rate.get_value(index));
          } else {
              float rate =  death_rate.get_value(index) / expected_max;
              row[local_index] = get_color(focal_cell_type, rate);
          }
      }
  }
}



void MainWindow::update_image(size_t sq_size,
                              const std::array< binned_distribution, 4 >& growth_rate) {

    cell_type focal_cell_type = normal;
    if( focal_display_type == normal_rate)    focal_cell_type = normal;
    if( focal_display_type == cancer_rate)    focal_cell_type = cancer;
    if( focal_display_type == infected_rate)  focal_cell_type = infected;
    if( focal_display_type == resistant_rate) focal_cell_type = resistant;

    if(grid_type == regular) {
      if(focal_display_type == dominant_rate) {
          display_regular(growth_rate);
      } else if (focal_display_type == cancer_death_rate) {
          display_regular_death_rate(growth_rate[cancer], cancer, sq_size);
      } else if (focal_display_type == normal_death_rate) {
          display_regular_death_rate(growth_rate[normal], normal, sq_size);
      } else {
          display_regular(growth_rate[focal_cell_type], focal_cell_type);
      }
    }

    if(grid_type == voronoi) {
        if(focal_display_type != dominant_rate) {
          display_voronoi(growth_rate[focal_cell_type], focal_cell_type, sq_size);
        } else {
          display_voronoi(growth_rate, sq_size);
        }
    }

    int w = ui->q_label->width();
    int h = ui->q_label->height();

    ui->q_label->setPixmap((QPixmap::fromImage(image_)).scaled(w,h, Qt::KeepAspectRatio));
    ui->q_label->update();
}

std::string get_string(std::string s, float v) {
    std::string output = s + " " + std::to_string(v) + "\n";
    return output;
}

void MainWindow::print_params(const Param& p) {
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
    s << get_string("Birth Rate Resistant", p.birth_cancer_resistant);
    s << get_string("Death Rate Resistant", p.death_cancer_resistant);

    s << get_string("Infection routine", p.infection_type);
    s << get_string("Infection %", p.percent_infected);
    s << get_string("Start type", p.start_setup);
    s << get_string("Prob normal infection", p.prob_normal_infection);
    s << get_string("Frequency resistant", p.freq_resistant);

    s << get_string("Distance infection upon death", p.distance_infection_upon_death);
    s << get_string("Probability infection upon death", p.prob_infection_upon_death);

    ui->text->appendPlainText(QString::fromStdString(s.str()));
    return;
}


void MainWindow::update_parameters(Param& p) {

   p.maximum_time = static_cast<int>(ui->box_maxtime->value());
   p.time_adding_cancer = static_cast<int>(ui->box_cancer_time->value());
   p.time_adding_virus = static_cast<int>(ui->box_virus_time->value());

   p.initial_number_cancer_cells = static_cast<size_t>(ui->box_cancer_cells->value());
   p.initial_number_normal_cells = static_cast<size_t>(ui->box_normal_cells->value());

   p.birth_normal = static_cast<float>(ui->box_birth_normal->value());
   p.death_normal = static_cast<float>(ui->box_death_normal->value());

   p.birth_cancer = static_cast<float>(ui->box_birth_cancer->value());
   p.death_cancer = static_cast<float>(ui->box_death_cancer->value());

   p.birth_infected = static_cast<float>(ui->box_birth_infected->value());
   p.death_infected = static_cast<float>(ui->box_death_infected->value());

   p.birth_cancer_resistant = static_cast<float>(ui->box_birth_cancer_resistant->value());
   p.death_cancer_resistant = static_cast<float>(ui->box_death_cancer_resistant->value());

   p.percent_infected = static_cast<float>(ui->box_percent_infected->value());
   p.prob_normal_infection = static_cast<float>(ui->box_prob_normal_infection->value());
   p.freq_resistant = static_cast<float>(ui->box_freq_resistant_cancer->value());

   p.distance_infection_upon_death = static_cast<float>(ui->box_distance_infection_death->value());
   p.prob_infection_upon_death = static_cast<float>(ui->box_prob_infection_death->value());

   p.diffusion = static_cast<float>(ui->box_diffusion->value());
   p.evaporation = static_cast<float>(ui->box_evaporation->value());
   p.t_cell_increase = static_cast<float>(ui->box_inflammation->value());
   p.t_cell_rate = static_cast<float>(ui->box_t_cell_rate->value());
   p.t_cell_density_scaler = static_cast<float>(ui->box_density_scaler->value());

   p.sq_num_cells = static_cast<size_t>(ui->box_sq_num_cells->value());

   p.infection_type = random_infection;

   auto infection_string = ui->box_infection_routine->currentText();
   if(infection_string == "Random")
       p.infection_type = random_infection;
   if(infection_string == "Center")
       p.infection_type = center_infection;
   if(infection_string == "Periphery")
       p.infection_type = periphery_infection;


    auto start_string = ui->box_start_setup->currentText();
    if(start_string == "Grow")
        p.start_setup = grow;
    if(start_string == "Full")
        p.start_setup = full;

   auto display_string = ui->drpdwnbox_display->currentText();
   if(display_string == "Cell types")
       focal_display_type = cells;
   if(display_string == "Normal Growth Rate")
       focal_display_type = normal_rate;
   if(display_string == "Normal Death Rate")
       focal_display_type = normal_death_rate;
   if(display_string == "Cancer Growth Rate")
       focal_display_type = cancer_rate;
   if(display_string == "Cancer Death Rate")
       focal_display_type = cancer_death_rate;
   if(display_string == "Infected Growth Rate")
       focal_display_type = infected_rate;
   if(display_string == "Resistant Growth Rate")
       focal_display_type = resistant_rate;
   if(display_string == "Dominant Growth Rate")
       focal_display_type = dominant_rate;
   if(display_string == "T-cells")
       focal_display_type = t_cells;

   auto grid_string = ui->box_grid_type->currentText();
   if(grid_string == "regular") {
       grid_type = regular;        // plotting flag
       p.use_voronoi_grid = false; // simulation flag
   }
   if(grid_string == "voronoi") {
       grid_type = voronoi;       // plotting flag
       p.use_voronoi_grid = true; // simulation flag
   }

   p.sq_num_pixels = static_cast<size_t>(ui->box_sq_num_pixels->value());

   if(p.use_voronoi_grid && p.sq_num_pixels <= p.sq_num_cells) {
       QMessageBox::warning(this,
                            tr("Oncolytic Virus Simulator"),
                            tr("Displaying Voronoi only works when\n"
                               "using more pixels than cells!"));
   }

   set_resolution(static_cast<int>(p.sq_num_pixels),
                  static_cast<int>(p.sq_num_pixels));

   print_params(p);
   return;
}

void MainWindow::update_polygons(const std::vector< std::vector< voronoi_point > >& all_edges) {

  factor_x = static_cast<float>(1.f * row_size / sim->sq_size);
  factor_y = static_cast<float>(1.f * col_size / sim->sq_size);

  polygons.clear();
  for(const auto& i : all_edges) {
      QPolygonF polygon;
      for(auto j : i) {
          polygon << QPointF(static_cast<qreal>(j.x_ * factor_x),
                             static_cast<qreal>(j.y_ * factor_y));
      }
      polygons.push_back(polygon);
  }
}

void MainWindow::setup_simulation() {
    if(is_running) {
        QMessageBox::warning(this,
                             tr("Oncolytic Virus Simulator"),
                             tr("Simulation is still running, please stop first"));
        return;
    }

    x_t.clear();
    y_n.clear();
    y_c.clear();
    y_i.clear();
    y_r.clear();

    update_parameters(all_parameters);

    sim = std::make_unique<simulation>(all_parameters);   //simulation(all_parameters);

    //Simulation.initialize_network();
    std::vector< std::vector< voronoi_point > > all_polys;
    sim->initialize_network(all_polys);

    if(all_parameters.use_voronoi_grid == true) {
      update_polygons(all_polys);
      set_resolution(static_cast<int>(all_parameters.sq_num_pixels),
                   static_cast<int>(all_parameters.sq_num_pixels));
    }
    if(all_parameters.use_voronoi_grid == false) {
        set_resolution(static_cast<int>(all_parameters.sq_num_cells),
                     static_cast<int>(all_parameters.sq_num_cells));
    }

    sim->t = 0.0;

    ui->btn_start->setText("Start");

    if(focal_display_type == cells) {
        update_image(all_parameters.sq_num_cells, false);
    } else if (focal_display_type == t_cells) {
        update_image(all_parameters.sq_num_cells, true);
    } else if (focal_display_type == cancer_death_rate) {
        update_image(all_parameters.sq_num_cells, sim->death_prob);
    } else if (focal_display_type == normal_death_rate) {
        update_image(all_parameters.sq_num_cells, sim->death_prob);
    } else {
        update_image(all_parameters.sq_num_cells, sim->growth_prob);
    }

    update_plot(static_cast<double>(sim->t),
                sim->get_count_cell_types());

    is_paused = true;

    QApplication::processEvents();
}



void MainWindow::on_btn_start_clicked()
{
    if(is_running) return; // pre emptive return if the simulation is already running

    if(!is_paused) setup_simulation(); // only setup if it is not paused

    is_paused = false; // now everything is setup, no need to pause.

    auto start_time = std::chrono::high_resolution_clock::now();
    ui->btn_start->setText("Start");

    is_running = true;
    int counter = 0;
    while(sim->t < all_parameters.maximum_time) {
        sim->update_one_step();
        counter++;

        int progress = static_cast<int>(100.f * sim->t / all_parameters.maximum_time);
        ui->progressBar->setValue(progress);

        // update speed is in 1 - 100
        const static size_t range = 1 * sim->world.size();

        int update_step = 1 + static_cast<int>((update_speed - 1) * 0.01f * range);

        if(counter % update_step == 0) {
            if(focal_display_type == cells) {
                update_image(all_parameters.sq_num_cells, false);
            } else if (focal_display_type == t_cells) {
                update_image(all_parameters.sq_num_cells, true);
            } else if (focal_display_type == cancer_death_rate) {
                update_image(all_parameters.sq_num_cells, sim->death_prob);
            } else if (focal_display_type == normal_death_rate) {
                update_image(all_parameters.sq_num_cells, sim->death_prob);
            } else {
                update_image(all_parameters.sq_num_cells, sim->growth_prob);
            }

            update_plot(static_cast<double>(sim->t),
                        sim->get_count_cell_types());
            QApplication::processEvents();
        }
        if(!is_running) break;
    }


    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( end_time - start_time ).count();

    float time_taken = 1.f * duration / 1000; // this is a rather ugly translation of milli seconds to seconds

    std::stringstream st;
    st << "This took ";
    st << time_taken;
    st << " seconds\n";
    ui->text->appendPlainText(QString::fromStdString(st.str()));
    std::cout << time_taken<< "\n";
    is_running = false;
}

void MainWindow::update_plot(double t,
                             const std::array<size_t, 5>& cell_numbers) {
    x_t.append(t);
    y_n.append(cell_numbers[0]);
    y_c.append(cell_numbers[1]);
    y_i.append(cell_numbers[2]);
    y_r.append(cell_numbers[3]);

    ui->line_plot->graph(0)->clearData();
    ui->line_plot->graph(0)->setData(x_t, y_n);

    ui->line_plot->graph(1)->clearData();
    ui->line_plot->graph(1)->setData(x_t, y_c);

    ui->line_plot->graph(2)->clearData();
    ui->line_plot->graph(2)->setData(x_t, y_i);

    ui->line_plot->graph(3)->clearData();
    ui->line_plot->graph(3)->setData(x_t, y_r);


    ui->line_plot->rescaleAxes();
    ui->line_plot->replot();
}

void MainWindow::on_btn_stop_clicked() {
    is_running = false;
    is_paused = true;
    ui->btn_start->setText("Resume");
}

void MainWindow::on_speed_slider_actionTriggered(int action) {
   update_speed = ui->speed_slider->value();
}

void MainWindow::on_drpdwnbox_display_activated(int index)
{
    auto display_string = ui->drpdwnbox_display->currentText();
    if(display_string == "Cell types")
        focal_display_type = cells;
    if(display_string == "Normal Growth Rate")
        focal_display_type = normal_rate;
    if(display_string == "Normal Death Rate")
        focal_display_type = normal_death_rate;
    if(display_string == "Cancer Growth Rate")
        focal_display_type = cancer_rate;
    if(display_string == "Cancer Death Rate")
        focal_display_type = cancer_death_rate;
    if(display_string == "Infected Growth Rate")
        focal_display_type = infected_rate;
    if(display_string == "Resistant Growth Rate")
        focal_display_type = resistant_rate;
    if(display_string == "Dominant Growth Rate")
        focal_display_type = dominant_rate;
    if(display_string == "T-cells")
        focal_display_type = t_cells;

    if (!is_running) {
        if(focal_display_type == cells) {
            update_image(all_parameters.sq_num_cells, false);
        } else if (focal_display_type == t_cells) {
            update_image(all_parameters.sq_num_cells, true);
        } else if (focal_display_type == cancer_death_rate) {
            update_image(all_parameters.sq_num_cells, sim->death_prob);
        } else if (focal_display_type == normal_death_rate) {
            update_image(all_parameters.sq_num_cells, sim->death_prob);
        } else {
            update_image(all_parameters.sq_num_cells, sim->growth_prob);
        }
          QApplication::processEvents();
    }
}

void MainWindow::on_btn_setup_clicked()
{
  ui->progressBar->setValue(0);
  setup_simulation();
}

void MainWindow::on_btn_add_virus_clicked()
{
   auto infection_string = ui->box_infection_routine->currentText();
   if(infection_string == "Random")
       sim->set_infection_type(random_infection);
   if(infection_string == "Center")
       sim->set_infection_type(center_infection);
   if(infection_string == "Periphery")
       sim->set_infection_type(periphery_infection);

   sim->set_percent_infected(static_cast<float>(ui->box_percent_infected->value()));

   sim->add_infected(sim->get_infection_type(),
                     sim->get_percent_infected());
}
