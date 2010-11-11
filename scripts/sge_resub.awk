BEGIN {
   option["sge_o_shell"] = "-S";
   option["sge_o_workdir"] = "-wd";
   option["stderr_path_list"] = "-e";
   option["stdout_path_list"] = "-o";
   option["hard resource_list"] = "-l";
   option["job_name"] = "-N";
   option["hard_queue_list"] = "-q";
   option["parallel environment"] = "-pe";
   option["jid_predecessor_list"] = "-hold_jid";
}

function parse_qstat(cmd) {

   delete job_info;
   delete error_task_ids;
   delete successors;
   delete succ_holds;

   while (cmd | getline line > 0) {
      split(line, field, ":  +");
      if (length(field) == 2) {
         if (field[1] ~ /^std..._path_list$/) {
            split(field[2], paths, ":");
            job_info[field[1]] = paths[3];
         }
         else if (field[1] ~ /^error reason/) {
            split(field[1], words);
            error_task_ids[words[3]] = field[2];
         }
         else if (field[1] == "jid_successor_list") {
            split(field[2], successors, ",");
            for (s in successors) {
               succ_cmd = "qstat -j " successors[s] " | grep \"jid_predecessor_list:\"";
               succ_cmd | getline sh;
               close(succ_cmd)
               split(sh, sh_fields, ":  +");
               succ_holds[successors[s]] = sh_fields[2];
            }
         }
         else if (field[1] == "parallel environment") {
            split(field[2], pe, ": ");
            job_info[field[1]] = "serial " pe[2];
         }
         else
            job_info[field[1]] = field[2];
      }
   }
   close(cmd)
}

{ 
   printf("processing %s            \r", $1) > "/dev/stderr"; 

   # parse qstat -j
   cmd = "qstat -j " $1; 
   parse_qstat(cmd);

   # build qsub cmd
   qsub_cmd = "qsub -terse -h";
   for (param in job_info) {
      if (option[param]) {
         qsub_cmd = qsub_cmd " " option[param] " " job_info[param];
      }
   }

   delete job_ids;

   # if this is a job array, resubmit one job per failed task id
   if (job_info["job-array tasks"]) {
      for (error_task in error_task_ids) {
         run_cmd = qsub_cmd " -t " error_task " " job_info["script_file"];
         run_cmd | getline job_id;
         close(run_cmd);
         print "Resubmitted " job_info["job_number"] "." error_task " as " job_id " (with hold)";
         job_ids[job_id] = 1;
      }
   }
   else {
      run_cmd = qsub_cmd " " job_info["script_file"];
      run_cmd | getline job_id;
      close(run_cmd);
      print "Resubmitted " job_info["job_number"] " as " job_id " (with hold)";
      job_ids[job_id] = 1;
   }

   for (successor in succ_holds) {
      jids = succ_holds[successor];
      for (job_id in job_ids) {
         jids = jids "," job_id;
      }
      qalter_cmd = "qalter " successor " -hold_jid " jids;
      qalter_cmd | getline qalter_status;
      close(qalter_cmd);
      print "Successor " successor " now depends on " jids ": " qalter_status;
   }

   if (job_info["job-array tasks"]) {
      for (error_task in error_task_ids) {
         qdel_cmd = "qdel " job_info["job_number"] " -t " error_task;
         qdel_cmd | getline qdel status;
         print "Deleted " job_info["job_number"] "." error_task ": " status;
      }
   }
   else {
      qdel_cmd = "qdel " job_info["job_number"];
      qdel_cmd | getline qdel status;
      print "Deleted " job_info["job_number"] ": " status;
   }

   for (job_id in job_ids) {
      qalter_cmd = "qalter " job_id " -h U";
      qalter_cmd | getline qalter_status;
      print "Released hold on " job_id ": " qalter_status;
   }

}
