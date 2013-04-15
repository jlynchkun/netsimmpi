%
% parentDirPath is usually
%   /share/data/bsc/fmri/sni/tal_mean100_coefvar0p05_head_g_t_148/corr-areas
% fmriDirPattern is usually one of the following
%   fmri_.*_ws
%   fmri_.*_ptsd
%   fmri_.*
% for testing specify the optional third argument -- for example 10
% for production do not specify the third argument
function s = allSubjectsAreaPairStatistics(parentDirPath, fmriDirPattern, optionalSubjectFileCountLimit)
  if nargin > 2
    subjectFileCountLimit = optionalSubjectFileCountLimit;
    display(['file(s) per subject limited to ' num2str(subjectFileCountLimit) ' file(s)'])
  else
    subjectFileCountLimit = Inf;
  end

  areaCount = 110;
  subjectDirPathList = java.util.LinkedList;
  % get a directory listing of the parent directory
  parentDirList = dir(parentDirPath);
  for subjectDirListItem = parentDirList'
    if subjectDirListItem.isdir 
      fmriDirPatternMatch = regexp(subjectDirListItem.name, fmriDirPattern, 'match');
      if isempty(fmriDirPatternMatch)
        display(['directory ' subjectDirListItem.name ' does not match directory pattern ' fmriDirPattern])
      else
        display(['found subject directory: ' subjectDirListItem.name])
        subjectDirPathList.add(fullfile(parentDirPath, subjectDirListItem.name));
      end
    else
      display(['found file: ' subjectDirListItem.name])
    end
  end
  
  subjectCount = subjectDirPathList.size();
  %s = zeros(areaCount, areaCount, subjectCount);
  % s(:,:,1) is for mean
  % s(:,:,2) is for standard deviation
  % s(:,:,3) is for sum(zcc0)
  % s(:,:,4) is for size(zcc0)
  s = nan(areaCount, areaCount, 2);
  display(['found ' num2str(subjectCount) ' subject directories matching directory pattern ' fmriDirPattern])

  %subjectIndex = 0;
  subjectDirPathListItr = subjectDirPathList.listIterator();
  while subjectDirPathListItr.hasNext()
    tic
    %subjectIndex = subjectIndex + 1;
    subjectDirPath = subjectDirPathListItr.next();
    subjectResultsDirPath = fullfile(subjectDirPath, 'results');
    display(['reading files in ' subjectDirPath])
    % read all corr-areas files in the subject directory
    subjectResultsDirList = dir(subjectResultsDirPath);
    corrAreasFileCount = 0;
    s(:) = nan;
    for subjectResultsDirListItem = subjectResultsDirList'
      if subjectResultsDirListItem.isdir
        % ignore directories
      elseif corrAreasFileCount >= subjectFileCountLimit
        display(['file limit ' num2str(subjectFileCountLimit) ' has been reached'])
        break
      else
        corrAreasFileCount = corrAreasFileCount + 1;
        corrAreasFilePath = fullfile(subjectResultsDirPath, subjectResultsDirListItem.name);
        c = bscReadCorrAreasFile(corrAreasFilePath);
        % sum zero-lag correlations
        zcc0 = atanh(c.voxelPairCorrelations(c.maxLag+1,:));
        cc0mean = tanh(mean(zcc0));
        cc0std = tanh(std(zcc0));
        zcc0sum = sum(zcc0);
        zcc0count = numel(zcc0);
        s(c.areaOne, c.areaTwo, 1) = cc0mean;
        s(c.areaTwo, c.areaOne, 1) = cc0mean;
        s(c.areaOne, c.areaTwo, 2) = cc0std;
        s(c.areaTwo, c.areaOne, 2) = cc0std;
        s(c.areaOne, c.areaTwo, 3) = zcc0sum;
        s(c.areaTwo, c.areaOne, 3) = zcc0sum;
        s(c.areaOne, c.areaTwo, 4) = zcc0count;
        s(c.areaTwo, c.areaOne, 4) = zcc0count;
      end
    end
    toc
    display(['read ' num2str(corrAreasFileCount) ' file(s) in ' subjectResultsDirPath])
    [~, subjectDirName] = fileparts(subjectDirPath);
    subjectCc0StatisticsFileName = ['cc0stats_' subjectDirName '.mat'];
    subjectCc0StatisticsFilePath = fullfile(subjectDirPath, subjectCc0StatisticsFileName);
    display(['writing ' subjectCc0StatisticsFilePath]);
    save(subjectCc0StatisticsFilePath, 's');
  end
end