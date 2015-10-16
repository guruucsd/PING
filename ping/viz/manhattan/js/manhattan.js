var functionOf = "";
var yvalue = "";
var analysis_names = [];
var version = "";
var masterlist = [];
var snps = [];

var ManhattanPlot = function(args) {
    var _this = this;
    _this.divID = args.divID;
    _this.gwas_url = args.gwas_url;
    _this.chrom_url = args.chrom_url || "data/chrom_locs.json";
    _this.snps_url = args.snps_url || "data/snps_all.json";
    _this.markerSize = args.markerSize || 4;

    _this.load_chrom_locations = function() {
        $.ajax({
            dataType: "json",
            url: _this.chrom_url,
            data: function(data) {},
            error: function(err) { alert('Load chrom error'); },
            success: function(data, textStatus, jqXHR) {
                _this.chromLocs = data;
                console.log('chromLocs');
                _this.load_snps()
            }
        })
    }

    _this.load_snps = function() {
        $.ajax({
            dataType: "json",
            url: _this.snps_url,
            data: function(data) {},
            error: function(err) { alert('Load snp error'); },
            success: function(data, textStatus, jqXHR) {
                _this.snps = data;
                console.log('snps');
                _this.load_gwas()
            }
        })

    _this.load_gwas = function() {
        $.ajax({
            dataType: "json",
            url: _this.gwas_url,
            data: function(data) {},
            error: function(err) { alert('Load error'); },
            success: function(data, textStatus, jqXHR) {
                _this.id = Object.keys(data)[0];
                _this.gwas = data[_this.id];
                console.log('gwas');
                _this.draw(_this.id, _this.chromLocs, _this.snps, _this.gwas);
            }
        });
    }

    _this.draw = function(id, chromLocs, snps, gwas) {
        // add a highchart to this id
        var thresh = 5e-8;

        // create x/y array
        var dataLabel = {
            data: [],
            showInLegend: 0,
            dataLabels: {
                enabled: true,
                x: 35,
                formatter: function() { return this.point.name; },
                style: { color: "white" }
            }
        };

        // Get max value, use parseFloat
        var minP = Object.keys(gwas).reduce(function (curVal, newVal) {
            return Math.min(curVal, gwas[newVal].pval);
        }, 1);

        var vals = Object.keys(gwas).map(function(curKey, idx) {
          var curVal = gwas[curKey];
          var chromIdx = parseInt(curVal.metadata.chromosome);
          var xval = chromLocs[chromIdx] + parseInt(curVal.metadata.basepair);
          var yval = -Math.log(curVal.pval) / Math.log(10);

          if (idx < 3) {
              // Assuming these are ordered (bad...),
              // show a laebel for it.
              dataLabel.data.push({
                  "x": xval,
                  "y": yval,
                  "z": 100,
                  "name": curKey
              });
          }
          return [xval, yval];
        })

        // now draw the vals array into #id
        var options = {
            chart: {
                renderTo: this.divID,
                defaultSeriesType: 'scatter',
                zoomType: 'x',
                backgroundColor: 'rgba(0,0,0,0)'
            },
            title: { enabled: true, text: 'title' },
            subtitle: { enabled: true, text: 'subtitle' },
            xAxis: {
                title: {
                    enabled: true,
                    text: 'chromosom'
                },
                tickWidth: 1,
                gridLineWidth: 0,
                labels: {
                    align: 'left',
                    x: 3,
                    y: 14
                }
            },
            yAxis: [{
                title: {
                    text: '-log<sub>10</sub> P',
                    useHTML: true
                },
                lineWidth: 1,
                tickWidth: 1,
                gridLineWidth: 0,
            }],
            legend: {
              align: 'right',
              verticalAlign: 'top',
              y: 45,
              floating: true,
              borderWidth: 1
            },
            tooltip: {
                formatter: function() {
                  var series = chart.get(this.series.name);
                  if (this.series.name == 'required') {
                    return 'By convention significance is measures as being above this line at -log<sub>10</sub> ' + thresh + '.';
                  }
                  if (this.series.name == 'floor') {
                    return 'Only SNPs with a value above -log<sub>10</sub> ' + minP + ' are displayed.';
                  }

                  // Find the matching point
                  var snp = null;
                  for (u = 0; u < this.series.data.length; u++) {
                      if (this.series.data[u].x == this.x && this.series.data[u].y == this.y) {
                          snp = Object.keys(gwas)[u];
                          break;
                        }
                  }

                  return snp + ' (chr.: ' + gwas[snp].metadata.chromosome + ')<br/>-log<sub>10</sub>P = ' + this.y.toFixed(2);
                },
                useHTML: true
            },
            plotOptions: {
              scatter: {
                point: {
                  events: {
                    touched: function() {
                      $(document).trigger('click', this);
                    },
                    click: function() {
                      console.log(this);
                    }
                  }
                },
                marker: {
                  radius: _this.markerSize,
                  states: {
                    hover: {
                      enabled: true,
                      lineColor: 'rgb(200,20,20)'
                    }
                  }
                },
                states: {
                  hover: {
                    marker: {
                      enabled: false
                    }
                  }
                }
              }
            }
        }; // end options

        // add the plot bands that signal the different chromosomes
        options.xAxis.plotBands = [];
        for (var i = 0; i < chromLocs.length - 1; i++) {
          options.xAxis.plotBands.push({
            from: chromLocs[i],
            to: chromLocs[i + 1],
            color: i % 2 == 0 ? 'rgba(255,127,80,.5)' : 'rgba(165,42,42,.5)',
            label: {
              text: (i == 0 || i == chromLocs.length - 2) ? '' : i,
              style: {
                color: '#606060'
              },
              verticalAlign: 'bottom',
              y: -9
            }
          });
        }

        var e = chromLocs[chromLocs.length - 1];
        // options.series[3] = dataLabel;

        options.series = [{          // What's this?
            data: vals,  //        options.series[0].data.dataLabels = dataLabel;
            name: 'genetics',
            id: 'genetics',
            marker: {
                radius: _this.markerSize,
                lineColor: 'rgba(10,10,100,1)',
                lineWidth: 1,
                fillColor: 'rgba(251,180,0,.6)'
            },
            shadow: {
                color: '#FF0000',
                offsetX: 0,
                offsetY: 0,
                width: 20
            },
            showInLegend: 0,
            zIndex: 1
        }, {  // lets add two lines, one for minimum displayed and one for required (red) at -log10(5e-8)
            data: [
              [0, vals[vals.length - 1][1]],
              [e, vals[vals.length - 1][1]]
            ],
            type: 'line',
            name: 'floor',
            id: 'floor',
            lineColor: 'rgba(100,100,100,.3)',
            lineWidth: .5,
            marker: {
                enabled: false
            },
            showInLegend: 0,
            zIndex: -1,
            enableMouseTracking: true
        }, {
          data: [
            [0, -Math.log(thresh) / Math.log(10)],
            [e, -Math.log(thresh) / Math.log(10)]
          ],
          type: 'line',
          name: 'required',
          id: 'required',
          lineColor: 'rgba(255,10,10,1)',
          lineWidth: 0.5,
          marker: {
            enabled: false
          },
          showInLegend: 0,
          zIndex: -1,
          enableMouseTracking: true
        }];
        options.series[3] = dataLabel;

        chart = new Highcharts.Chart(options);
    } // end .draw
}
    _this.load_chrom_locations()
}