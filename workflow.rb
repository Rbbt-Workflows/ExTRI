require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

require 'rbbt/sources/ExTRI'
require 'rbbt/sources/HTRI'
require 'rbbt/sources/TRRUST'
require 'rbbt/sources/tfacts'
require 'rbbt/sources/TFCheckpoint'
require 'rbbt/sources/signor'
require 'rbbt/sources/GOA'
require 'rbbt/sources/Intact'
require 'rbbt/sources/uniprot'
require 'rbbt/sources/CytReg'
require 'rbbt/sources/GEREDB'
require 'rbbt/sources/Pavlidis'

Workflow.require_workflow "Appris"
module ExTRI
  extend Workflow
end

require 'rbbt/tasks/ExTRI/clean'
require 'rbbt/tasks/ExTRI/validation'
#require 'rbbt/tasks/ExTRI/validation_samples'
require 'rbbt/tasks/ExTRI/validation_evaluation'
require 'rbbt/tasks/ExTRI/database_coverage'
require 'rbbt/tasks/ExTRI/ner'
require 'rbbt/tasks/ExTRI/statistics'
require 'rbbt/tasks/ExTRI/year'
require 'rbbt/tasks/ExTRI/psicquic'
require 'rbbt/tasks/ExTRI/biogateway'
require 'rbbt/tasks/ExTRI/regulon'
require 'rbbt/tasks/ExTRI/greco'
require 'rbbt/tasks/ExTRI/tf_tf'
require 'rbbt/tasks/ExTRI/plots'
require 'rbbt/tasks/ExTRI/adhoc'
require 'rbbt/tasks/ExTRI/manuscript'
require 'rbbt/tasks/ExTRI/rename'

ExTRI.export :CollecTRI, :ExTRI_final

