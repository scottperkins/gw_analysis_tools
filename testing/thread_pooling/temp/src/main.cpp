#include <iostream>
//#include "test.h"
//#include "test2.h"
//#include "wrapper.h"
//#include "wrapper2.h"
#include <adolc/adouble.h>
#include <cstdlib>
#include <string>
#include <vector>
#include <ctime>
#include <numeric>
#include <cmath>
#include <sstream>
#include <thread>
#include <chrono>
#include <mutex>
#include <unistd.h>
#include <functional>
#include <condition_variable>
#include <queue>
#include <omp.h>
#include <time.h>

using namespace std;

struct sampler
{
	double *data;
	double *results;
	bool *finished;
	int data_length;	
	int *counter;

};


sampler *samplerptr;


void do_something(sampler *sampler, int id)
{
	for (int j =0 ; j<(10); j++)
	{
	double result=0;
	for (int i =0; i< sampler->data_length; i++)
		result += log(sampler->data[i]/1.e3); 
	sampler->results[id] = result;
	}
	sampler->finished[id] = true;
}

void do_something3(int id)
{
	//std::cout<<"DOING SOMETHING"<<std::endl;
	do_something(samplerptr, id);
	samplerptr->counter[id]+=1;
	samplerptr->finished[id] = true;
}

void do_something2(int id)
{
	std::cout<<id<<std::endl;
}

class ThreadPool

{
public:
	using Task = std::function<void()>;
	
	
	explicit ThreadPool(std::size_t numThreads)
	{
		start(numThreads);
	}

	~ThreadPool()
	{
		stop();
	}


	//void enqueue(Task task)
	//{
	//	{
	//		std::unique_lock<std::mutex> lock{mEventMutex};
	//		mTasks.emplace(std::move(task));
	//	}
	//	mEventVar.notify_one();
	//}
	void enqueue(int i)
	{
		{
			std::unique_lock<std::mutex> lock{mEventMutex};
			mTasks.emplace(std::move(i));
		}
		mEventVar.notify_one();
	}
	
	void public_stop()
	{
		stop();
	}
private:
	std::vector<std::thread> mThreads;
	
	std::condition_variable mEventVar;

	std::mutex mEventMutex;

	bool mStopping = false;
		
	//std::queue<Task> mTasks;
	std::queue<int> mTasks;

	void start(std::size_t numThreads)
	{
		for(auto i =0u; i<numThreads; i++)
		{
			mThreads.emplace_back([=]{
				while(true)
				{
					int j;
					{
						std::unique_lock<std::mutex> lock{mEventMutex};

						//mEventVar.wait(lock,[=]{return mStopping || !mTasks.empty(); });
						mEventVar.wait(lock,[=]{return mStopping || !mTasks.empty(); });
						
						if (mStopping && mTasks.empty())
							break;	
						j = std::move(mTasks.front());
						mTasks.pop();
						//std::cout<<mTasks.empty();
					}
					do_something3(j);
					
				}
			});
		}
	}
	void stop() noexcept
	{
		std::cout<<"STOPPING"<<std::endl;
		{
			std::unique_lock<std::mutex> lock{mEventMutex};
			mStopping = true;
		}
		
		mEventVar.notify_all();
		
		for(auto &thread: mThreads)
			thread.join();
	}
};

int main (void)
{
	int length_data = 1000;
	int chains= 100;
	int iterations = 1000;
	sampler samplerinst;
	samplerinst.data_length = length_data;
	samplerinst.data = (double *)malloc(sizeof(double)*length_data);
	for(int i =0 ; i<length_data;i++)
		samplerinst.data[i] = i;	
	samplerinst.results= (double *)malloc(sizeof(double)*chains);
	samplerinst.finished = (bool *)malloc(sizeof(bool)*chains);
	samplerinst.counter = (int *)malloc(sizeof(int)*chains);
	for(int i =0 ; i<chains;i++)
	{
		samplerinst.finished[i] = true;	
		samplerinst.counter[i]=0;
	}

	samplerptr = &samplerinst;
	
	int k=0, l =0;
	{
		ThreadPool pool(20);
		while(k < iterations && l<10*chains*iterations)
		{
			for (int i = 0; i<chains; i++){	
				if(samplerptr-> finished[i] && samplerptr->counter[i]<iterations){
					samplerptr->finished[i]=false;
					pool.enqueue(i);
					if(i ==0){k++;}//std::cout<<"K: "<<k<<std::endl;}
					l+=1;
					//std::cout<<l<<std::endl;
				}

			}
			//if(k>iterations || l>chains*iterations) pool.public_stop();
		}
		//pool.public_stop();
	}

		//ThreadPool pool(20);


		//pool.enqueue([] {
		//	std::cout<<"1"<<std::endl;
		//});

		//pool.enqueue([] {
		//	std::cout<<"2"<<std::endl;
		//});


	//while(k < iterations && l<chains*iterations)
	//{
	//	for (int i = 0; i<chains; i++){	
	//		if(samplerptr-> finished[i]){
	//			//pool.enqueue(i);
	//			do_something3(i);
	//			//samplerptr->finished[i]=false;
	//			if(i ==0){k++;}//std::cout<<"K: "<<k<<std::endl;}
	//			l+=1;
	//		}

	//	}
	//}

	for(int i =0; i< chains; i++)
	{
		std::cout<<samplerinst.results[i]<<std::endl;
		std::cout<<samplerinst.finished[i]<<std::endl;
		std::cout<<"COUNTER: "<<samplerinst.counter[i]<<std::endl;
	}
	free(samplerinst.data);
	free(samplerinst.results);
	free(samplerinst.finished);

	////####################################################

	//std::thread threads[chains];
	//int k=0;
	//while(k<iterations)
	//{
	//	for(int i =0; i< chains; i++)
	//	{
	//		//threads[i]= std::thread(do_something2, i);
	//		if(sampler.finished[i]){
	//			std::cout<<i<<std::endl;
	//			threads[i]= std::thread(do_something,&sampler, i);
	//			sampler.finished[i]=false;
	//			//threads[i].detach();
	//			if(i==0) {k++;std::cout<<"COUNTER: "<<k<<std::endl;}
	//		}
	//	}
	//	usleep(5e6);
	//	
	//}

	//for(int i =0; i< chains; i++)
	//	threads[i].join();

	////####################################################



	return 0;
}
