
classdef TestBenchmark_Agniva < matlab.unittest.TestCase
    
        properties (Access = public)
        
            data = '';
            test_rs_rectification_1 = '';
            
        end
        methods(TestMethodSetup)
            % Clearing all variables and closing all figures
            function clearVars(testCase)
                clear all;
                close all;
            end
            
            % Setting up data
            function setData(testCase)
                testCase.data = 'D:\code\CVPR_code_RS\CVPR_code\benchmark\SampleData\Yizhen_House_TestData'; % Change the target dataset here
                testCase.test_rs_rectification_1 = 'D:\code\CVPR_code_RS\CVPR_code\benchmark\SampleData\test_output.png';
            end
        end
    
    methods(TestMethodTeardown)
        % Clearing all variables and closing all figures again
        function clearVarsAndFigs(testCase)
            clear all;
            close all;
        end
    end
    
    methods(Test)
        
        % Basic Sanity check: that the function runs, does not crashes and
        % returns a 'TRUE' at the end of the program execution
        function basicSanity(testCase)
            addpath('./benchmark/');
            data = load(testCase.data);
            rs = BenchmarkYizhen;
            [rectifiedImg, poses,landmarkHomography, retVal] = rs.SolveImageRectification(data.IP1, data.IP2, data.K, data.K, data.image1, data.image2);
            
           verifyTrue(testCase,isequal(retVal,true));
        end
        
        % Structural Similarity Check: between two images, the first one
        % obtained from Yizhen's original code and test case, the second
        % one from our wrapper
        function structuralSimilarity(testCase)
            addpath('./benchmark/');
            data = load(testCase.data);
            rs = BenchmarkYizhen;
            [rectifiedImg, poses,landmarkHomography, retVal] = rs.SolveImageRectification(data.IP1, data.IP2, data.K, data.K, data.image1, data.image2);
            
            imgCompare = imread(testCase.test_rs_rectification_1);
            
            I1 = rgb2gray(rectifiedImg);
            I2 = rgb2gray(imgCompare);
            
            ssimVal = 1.0;
            
            fprintf('Obtained size of rectification image: (%d,%d); expected size: (%d,%d)\n', size(rectifiedImg,1), size(rectifiedImg,2), size(imgCompare,1), size(imgCompare,2));
            fprintf('Obtained aspect ratio of rectification image: (%d); expected size: (%d)\n\n', size(rectifiedImg,1)/size(rectifiedImg,2), size(imgCompare,1)/size(imgCompare,2));
            
            try
                [ssimVal,ssimMap] = ssim(double(I1),double(I2));
                fprintf('----> Obtained SSIM: %d\n', ssimVal);
            catch ME
                warning('Failed while computing SSIM, hence: RE-SIZING');
                I1 = imresize(I1, size(I2));
                [ssimVal,ssimMap] = ssim(double(I1),double(I2));
            end            
            
            figure;
            subplot(1,3,1); imshow(I1);
            subplot(1,3,2); imshow(I2);
            subplot(1,3,3); imshow(ssimMap,[]);
            title("Similarity index:"+num2str(ssimVal));
            
           verifyLessThan(testCase,ssimVal,0.2);
           
           drawnow();
           pause(7);
        end
        
        % Landmark Pose Check: simple sanity check for extracted landmark
        % poses
        function landmarkpose(testCase)
            addpath('./benchmark/');
            data = load(testCase.data);
            rs = BenchmarkYizhen;
            [rectifiedImg, poses, landmarkHomography, retVal] = rs.SolveImageRectification(data.IP1, data.IP2, data.K, data.K, data.image1, data.image2);
            
            pcshow(poses,'r','MarkerSize',100); hold on;
            xlabel('X'); ylabel('Y'); zlabel('Z');
            grid on;
            drawnow();
            pause(7);
            close all;            
            
           hasSomeNans = sum(sum(isnan(poses)));
           verifyEqual(testCase,hasSomeNans,0);

        end
        
        
    end
end