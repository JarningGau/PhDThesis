MeasureMultiRocAuc = R6::R6Class(
  "MeasureMultiRocAuc",
  inherit = mlr3::MeasureClassif,
  public = list(
    initialize = function() {
      super$initialize(
        # custom id for the measure
        id = "mclassif.mauc_aunu",

        # additional packages required to calculate this measure
        packages = character(),

        # properties, see below
        properties = character(),

        # required predict type of the learner
        predict_type = "prob",

        # feasible range of values
        range = c(0, 1),

        # minimize during tuning?
        minimize = TRUE
      )
    }
  ),

  private = list(
    # custom scoring function operating on the prediction object
    .score = function(prediction, ...) {
      mlr3measures::mauc_aunu(prediction$truth, prediction$prob)
    }
  )
)


mlr3::mlr_measures$add("mclassif.mauc_aunu", MeasureMultiRocAuc)

