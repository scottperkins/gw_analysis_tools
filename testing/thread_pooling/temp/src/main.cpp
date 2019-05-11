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
	int *swpcounter;

};


sampler *samplerptr;


void do_something(sampler *sampler, int id);
void do_something3(int id);
void do_something_else(int id);
void do_swap(int i, int j);
void do_something2(int id);
void do_something3_deterministic(int id);

class ThreadPool

{
public:
	//using Task = std::function<void()>;
	
	
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

	void enqueue_swap(int i)
	{
		{
			std::unique_lock<std::mutex> lock{mEventMutexSWP};
			mSwaps.emplace(std::move(i));
		}
		mEventVarSWP.notify_one();
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
	
	int numSwpThreads= 1;
		
	std::condition_variable mEventVarSWP;

	std::mutex mEventMutexSWP;

	//std::queue<Task> mTasks;
	std::queue<int> mTasks;
	std::queue<int> mSwaps;

	void start(std::size_t numThreads)
	{
		for(auto i =0u; i<numThreads-numSwpThreads; i++)
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

		//Swapping thread
		for(auto i =0u; i<numSwpThreads; i++)
		{
			mThreads.emplace_back([=]{
				while(true)
				{
					int j, k;
					{
						std::unique_lock<std::mutex> lock{mEventMutexSWP};

						//mEventVarSWP.wait(lock,[=]{return mStopping || !mSwaps.empty(); });
						mEventVarSWP.wait(lock,[=]{return mStopping || !(mSwaps.size()<2); });
						
						//if (mStopping && mSwaps.empty())
						if (mStopping && mSwaps.size()<2)
							break;	
						j = std::move(mSwaps.front());
						mSwaps.pop();
						k = std::move(mSwaps.front());
						mSwaps.pop();
						//std::cout<<mTasks.empty();
					}
					do_swap(j,k);
					
				}
			});
		}
	}
	void stop() noexcept
	{
		std::cout<<"STOPPING"<<std::endl;
		{
			std::unique_lock<std::mutex> lock{mEventMutex};
			//std::unique_lock<std::mutex> lock{mEventMutexSWP};
			mStopping = true;
		}
		
		mEventVar.notify_all();
		mEventVarSWP.notify_all();
		
		for(auto &thread: mThreads)
			thread.join();
	}
};
ThreadPool *poolptr;
int main (int argc, char *argv[])
{
	int option = atoi(argv[1]);

	int length_data = 1000;
	int chains= 20;
	int iterations = 10000;
	sampler samplerinst;
	samplerinst.data_length = length_data;
	samplerinst.data = (double *)malloc(sizeof(double)*length_data);
	for(int i =0 ; i<length_data;i++)
		samplerinst.data[i] = 1.+i;	
	samplerinst.results= (double *)malloc(sizeof(double)*chains);
	samplerinst.finished = (bool *)malloc(sizeof(bool)*chains);
	samplerinst.counter = (int *)malloc(sizeof(int)*chains);
	samplerinst.swpcounter = (int *)malloc(sizeof(int)*chains);
	for(int i =0 ; i<chains;i++)
	{
		samplerinst.finished[i] = true;	
		samplerinst.counter[i]=0;
		samplerinst.swpcounter[i]=0;
	}

	samplerptr = &samplerinst;
	
	int k=0, l =0;
	if(option==1)
	{
		ThreadPool pool(20);
		poolptr = &pool;	
		while(k < iterations && l<10*chains*iterations)
		{
			for (int i = 0; i<chains; i++){	
				//if(samplerptr-> finished[i] && samplerptr->counter[i]<iterations){
				if(samplerptr-> finished[i] ){
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



	if(option ==2)
	{
		while(k < iterations && l<chains*iterations)
		{
			for (int i = 0; i<chains; i++){	
				if(samplerptr-> finished[i]){
					//pool.enqueue(i);
					do_something3_deterministic(i);
					//samplerptr->finished[i]=false;
					//if(i ==0){k++;}//std::cout<<"K: "<<k<<std::endl;}
					//l+=1;
				}

			}
			for (int i = 0; i<chains; i++){	
				if(samplerptr-> finished[i]){
					//pool.enqueue(i);
					do_something_else(i);
					//samplerptr->finished[i]=false;
					if(i ==0){k++;}//std::cout<<"K: "<<k<<std::endl;}
					l+=1;
				}

			}
		}
	}

	for(int i =0; i< chains; i++)
	{
		std::cout<<"CHAIN: "<<i<<" | RESULT: "<<samplerinst.results[i]
			<< " | COUNTER: "<<samplerinst.counter[i]<< 
			" | SWPCOUNTER: "<<samplerinst.swpcounter[i]<<std::endl;
	}
	free(samplerinst.data);
	free(samplerinst.results);
	free(samplerinst.finished);

	//ThreadPool pool(20);


	//pool.enqueue([] {
	//	std::cout<<"1"<<std::endl;
	//});

	//pool.enqueue([] {
	//	std::cout<<"2"<<std::endl;
	//});

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
void do_something(sampler *sampler, int id)
{
	for (int j =0 ; j<(10); j++)
	{
	double result=0;
	for (int i =0; i< sampler->data_length; i++)
		result += log(sampler->data[i]*(id+1)/1.e3); 
	sampler->results[id] = result;
	}
	//sampler->finished[id] = true;
}

void do_something3(int id)
{
	//std::cout<<"DOING SOMETHING"<<std::endl;
	do_something(samplerptr, id);
	samplerptr->counter[id]+=1;
	poolptr->enqueue_swap(id);	
	//samplerptr->finished[id] = true;
}
void do_something3_deterministic(int id)
{
	//std::cout<<"DOING SOMETHING"<<std::endl;
	do_something(samplerptr, id);
	samplerptr->counter[id]+=1;
	//samplerptr->finished[id] = true;
}

void do_swap(int i, int j)
{
	//std::cout<<"DO SOMETHING ELSE THREAD "<<id<<std::endl;
	double temp = samplerptr->results[i];
	samplerptr->results[i] =samplerptr->results[j];	
	samplerptr->results[j] =temp;	
	samplerptr->swpcounter[i] +=1;
	samplerptr->swpcounter[j] +=1;
	samplerptr->finished[i] = true;
	samplerptr->finished[j] = true;
}
void do_something_else(int id)
{
	//std::cout<<"DO SOMETHING ELSE THREAD "<<id<<std::endl;
	samplerptr->results[id] = id;	
	samplerptr->finished[id] = true;
	samplerptr->swpcounter[id] +=1;
}
void do_something2(int id)
{
	std::cout<<id<<std::endl;
}
