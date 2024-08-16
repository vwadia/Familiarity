function manual_face_landmark_edit_gui

Initial_path = 'D:\Project2 memory\2D_human_face_model\faces in movies\images\*.*';
Initial_path = 'D:\Project2 memory\2D_human_face_model\familiar faces\images\*.*';

% space: load image
% return: save mark
% w,a,s,d: 
% right mouse: move
%
%


load landmarks_prototype.mat marks_ori

h = figure;
set(h, 'WindowButtonDownFcn',@mousebd)
set(h, 'WindowKeyPressFcn',@keydown)

im = [];
marks = [];
marks_old = [];
h_mark = [];
h_text = [];

xinit = [];
yinit = [];

mark_ind = [];
imagefilename = [];
isnew = 0;

imagefiles = find_files(Initial_path);
file_ind = 1;

load_image_and_mark()

    function mousebd(src,~)
        % mouse button down function
        % edit landmark
        
        if isempty(marks)
            return
        end
        
        cp = get(gca, 'CurrentPoint');
        xinit = cp(1,1);
        yinit = cp(1,2);
        marks_old = marks;
        
        seltype = src.SelectionType; %fprintf('%s\n',seltype)
        if strcmp(seltype,'normal')
            [~, mark_ind] = min( (xinit-marks(:,1)).^2 + (yinit-marks(:,2)).^2 );
        end
        
        src.WindowButtonMotionFcn = @mousebm;
        src.WindowButtonUpFcn = @mousebu;
    end

    function mousebm(src,~)
        cp = get(gca, 'CurrentPoint');
        x = cp(1,1);
        y = cp(1,2);

        seltype = src.SelectionType;
        if strcmp(seltype,'alt')
           marks(:,1) = marks_old(:,1) + x - xinit;
           marks(:,2) = marks_old(:,2) + y - yinit;
        elseif strcmp(seltype,'normal')
            marks(mark_ind,1) = x;
            marks(mark_ind,2) = y;
        end
        
        set(h_mark,'YData', marks(:,2))
        set(h_mark,'XData', marks(:,1))
        
    end

    function mousebu(src,~)
        src.WindowButtonMotionFcn = '';
        src.WindowButtonUpFcn = '';
    end


    function keydown(src,evnt)
        % fprintf('%s\n',evnt.Key)
        if strcmp(evnt.Key, 'space')
            
            % load image
            [FileName, PathName ] = uigetfile(Initial_path);
            
            imagefilename = fullfile(PathName,FileName);
            load_image_and_mark(imagefilename)
            
        elseif strcmp(evnt.Key, 'rightarrow')
            if isnew
                %save_mark_file
            end
            if file_ind< length(imagefiles)
                file_ind = file_ind+1;
                load_image_and_mark
            end
            
        elseif strcmp(evnt.Key, 'leftarrow')
            if isnew
                %save_mark_file
            end
            if file_ind > 1
                file_ind = file_ind-1;
                load_image_and_mark
            end
            
        elseif strcmp(evnt.Key, 'w')
            marks(:,2) = mark_scale(marks(:,2), 1.1);
            set(h_mark,'YData', marks(:,2))
            
        elseif strcmp(evnt.Key, 's')
            marks(:,2) = mark_scale(marks(:,2), 0.9);
            set(h_mark,'YData', marks(:,2))
            
        elseif strcmp(evnt.Key, 'a')
            marks(:,1) = mark_scale(marks(:,1), 0.9);
            set(h_mark,'XData', marks(:,1))
            
        elseif strcmp(evnt.Key, 'd')
            marks(:,1) = mark_scale(marks(:,1), 1.1);
            set(h_mark,'XData', marks(:,1))
            
        elseif strcmp(evnt.Key, 'z')
            marks = mark_rotate(marks, 1);
            set(h_mark,'XData', marks(:,1), 'YData', marks(:,2))
            
        elseif strcmp(evnt.Key, 'x')
            marks = mark_rotate(marks, -1);
            set(h_mark,'XData', marks(:,1), 'YData', marks(:,2))
            
        elseif strcmp(evnt.Key, 'return')
            save_mark_file
            msgbox('saved')
        end
        
    end


    function load_image_and_mark(filename)
        
        % load current image file
        if nargin<1
            imagefilename = imagefiles{file_ind};
        else
            imagefilename = filename;
        end
        
        if ~exist(imagefilename,'file')
            imagefilename = [];
            return
        end
        
%         if strcmp( FileName(end-2:end), 'mat')
%             % load other marks file for loaded image
%             matfilename = fullfile(PathName,FileName);
%             load(matfilename, 'marks');
%         end
        
        im = imread(imagefilename);
        cla
        image(im);
        axis equal
        hold on
        
        data_file = imagefile2markfile(imagefilename);
        if exist(data_file,'file') 
            % already labeled
            load(data_file,'pts'); marks = pts;
            isnew = 0;
            h_mark = plot(marks(:,1),marks(:,2),'r+-');
        else
            %marks = load_raw_mark_file( imagefilename );
            marks = marks_ori;
            isnew = 1;
            h_mark = plot(marks(:,1),marks(:,2),'g+-');
        end
        
%         marks = marks_ori;

                % label each point
%                 marks_no_nan = marks(~isnan(marks(:,1)),:); 
%                  
%                 for i = 1:length(marks_no_nan)
%                     mark_name{i} = num2str(i);
%                 end
%                 text(marks_no_nan(:,1),marks_no_nan(:,2),mark_name, 'color','w')
                
    end

    function save_mark_file()
        data_file = imagefile2markfile(imagefilename);
        pts = marks;
        save(data_file, 'pts');
    end

    function markfile = imagefile2markfile(imagefile)
        [path_str, name, ext] = fileparts(imagefile);
        markfile = fullfile(path_str, '..', 'landmarks80',  [name '.mat']);
    end


end


function x = mark_rotate(x, a)
a = a/180*pi;
m = x(~isnan(x(:,1)),:); %remove nan
m = m * [cos(a) -sin(a); sin(a) cos(a)] ;
x(~isnan(x)) = m;
end

function x = mark_scale(x, s)
m = nanmean(x);
x = x - m;
x = x * s+m;
end


function hm = draw_marks(marks)

end

