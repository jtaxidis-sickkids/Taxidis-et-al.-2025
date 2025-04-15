function [prediction,accuracy,svm] = SVM(Rtrain,Rpred,traintr,predtr)

svm = fitcsvm(Rtrain,traintr,'KernelFunction','rbf','Standardize',true,...
    'KernelScale','auto','ClassNames',[-1,1]);                              % TRAIN SVM MODEL

prediction = predict(svm,Rpred);                                            % PREDICT THE TRIALS
corrects = (predtr - prediction == 0);                                      % Find correct predictions
accuracy = 100*sum(corrects)/length(predtr);                                % Compute their percentage (Accuracy)