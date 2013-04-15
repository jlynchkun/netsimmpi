%
% parentDirPath is usually
%   ~/dev-local/netsimmpi/jobs/
% for testing specify the optional third argument -- for example 10
% for production do not specify the third argument
function s = allJobsStatistics(parentDirPath, jobDirPattern, optionalJobCountLimit)
  if nargin > 2
    jobCountLimit = optionalJobCountLimit;
    display(['jobs limited to ' num2str(jobCountLimit) ' jobs(s)'])
  else
    jobCountLimit = Inf;
  end

  areaCount = 110;
  jobDirPathList = java.util.LinkedList;
  % get a directory listing of the parent directory
  parentDirList = dir(parentDirPath);
  for jobDirListItem = parentDirList'
    if jobDirListItem.isdir 
      jobDirPatternMatch = regexp(jobDirListItem.name, jobDirPattern, 'match');
      if isempty(jobDirPatternMatch)
        display(['directory ' jobDirListItem.name ' does not match directory pattern ' jobDirPattern])
      else
        display(['found job directory: ' jobDirListItem.name])
        jobDirPathList.add(fullfile(parentDirPath, jobDirListItem.name));
      end
    else
      display(['found file: ' jobDirListItem.name])
    end
  end

  jobCount = jobDirPathList.size();
  %s = zeros(areaCount, areaCount, subjectCount);
  % s(:,:,1) is for mean
  % s(:,:,2) is for standard deviation
  % s(:,:,3) is for sum(zcc0)
  % s(:,:,4) is for size(zcc0)
  s = nan(areaCount, areaCount, 2);
  display(['found ' num2str(jobCount) ' directories matching directory pattern ' jobDirPattern])

  jobDirPathListItr = jobDirPathList.listIterator();
  while jobDirPathListItr.hasNext()
    tic
    jobDirPath = jobDirPathListItr.next();
    jobResultsDirPath = fullfile(jobDirPath, 'results');
    display(['reading files in ' jobDirPath])

    jobResultsDirList = dir(jobResultsDirPath);
    jobFileCount = 0;
    s(:) = nan;
    for jobResultsDirListItem = jobResultsDirList'
      if jobResultsDirListItem.isdir
        % ignore directories
      elseif jobFileCount >= jobCountLimit
        display(['file limit ' num2str(jobCountLimit) ' has been reached'])
        break
      else
        jobFileCount = jobFileCount + 1;
        resultsFilePath = fullfile(jobResultsDirPath, jobResultsDirListItem.name);
        %c = bscReadCorrAreasFile(resultsFilePath);
      end
    end
    toc
    display(['read ' num2str(jobFileCount) ' file(s) in ' jobResultsDirPath])
    [~, jobDirName] = fileparts(jobDirPath);
    jobStatisticsFileName = ['stats_' jobDirName '.mat'];
    jobStatisticsFilePath = fullfile(jobDirPath, jobStatisticsFileName);
    display(['writing ' jobStatisticsFilePath]);
    save(jobStatisticsFilePath, 's');
  end
end